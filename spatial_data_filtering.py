"""
Main RETA pipeline: QGIS Processing Algorithm and DataCleaner (V8 sequential).

Use this module when:
  - Running inside QGIS (Processing Toolbox) for full spatial data filtering.
  - You need the full pipeline: operational -> global -> spatial phases, smart turns,
    grid overlap, auto speed/swath/delay, and apply_filters_sequential().

For headless or lightweight use without QGIS, see data_cleaner.py (standalone
DataCleaner with apply_filters() and medcouple/AOI only).
"""

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterFeatureSink,
    QgsFeatureSink,
    QgsFeature,
    QgsField,
    QgsGeometry,
    QgsPointXY,
    QgsProcessingParameterDefinition
)
from qgis.PyQt.QtCore import QVariant
import pandas as pd
import numpy as np
import math

try:
    from .crs_utils import (
        is_geographic,
        distances_m,
        points_to_meters,
        area_m2_from_points,
        _normalize_crs,
    )
    CRS_UTILS_AVAILABLE = True
except ImportError:
    CRS_UTILS_AVAILABLE = False

# --- V8 SEQUENTIAL LOGIC (Robust Transects + Grid Overlap + Strict Phases) ---
class DataCleaner:
    def __init__(self, dataframe, crs=None, feedback=None):
        """
        dataframe: must have 'lat', 'lon' (and optionally 'time', variable column).
        crs: QgsCoordinateReferenceSystem or authid string (e.g. EPSG:4326). If set,
             distance/area use proper projection (no preset scale). If None, fallback
             uses WGS84 Haversine / approximate degree-to-m scaling.
        feedback: optional QgsProcessingFeedback for progress and cancel in long loops.
        """
        self.df = dataframe.copy()
        self._crs = (_normalize_crs(crs) if crs is not None else None) if CRS_UTILS_AVAILABLE else None
        self._feedback = feedback
        self._canceled = False

        # Taxonomy
        self.error_masks = {
            'Clean': 0,
            # Operational
            'Op: Turn': 1,
            'Op: Short Segment': 2,
            'Op: Distance Jump': 4,
            'Op: Collocated': 8,
            'Op: Start Delay': 16,
            'Op: End Delay': 32,
            'Op: Overlap': 64,
            # Global
            'Gl: Value Low': 128,
            'Gl: Value High': 256,
            # Local
            'Lo: Spatial Outlier': 512,
            # Limits
            'Op: Min Speed': 1024,
            'Op: Max Speed': 2048,
            'Lo: Cross-Track Outlier': 4096,
            'Op: High Acceleration': 8192,
            'Op: Variable Speed Change': 16384 # Inspired by YE Smooth Velocity
        }

    def calculate_geometry(self):
        self.df['prev_lat'] = self.df['lat'].shift(1)
        self.df['prev_lon'] = self.df['lon'].shift(1)

        if self._crs and CRS_UTILS_AVAILABLE:
            self.df['dist_m'] = distances_m(
                self.df['prev_lon'].values,
                self.df['prev_lat'].values,
                self.df['lon'].values,
                self.df['lat'].values,
                self._crs,
            )
            if is_geographic(self._crs):
                # Bearing from geographic: use standard formula
                dlambda = np.radians(self.df['lon'] - self.df['prev_lon'])
                phi1 = np.radians(self.df['prev_lat'])
                phi2 = np.radians(self.df['lat'])
                y = np.sin(dlambda) * np.cos(phi2)
                x = np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(dlambda)
                self.df['bearing_deg'] = np.degrees(np.arctan2(y, x))
            else:
                # Projected: bearing from dx, dy (e.g. UTM)
                dx = self.df['lon'] - self.df['prev_lon']
                dy = self.df['lat'] - self.df['prev_lat']
                self.df['bearing_deg'] = (90 - np.degrees(np.arctan2(dy, dx))) % 360
        else:
            R = 6371000.0
            phi1 = np.radians(self.df['prev_lat'])
            phi2 = np.radians(self.df['lat'])
            dphi = np.radians(self.df['lat'] - self.df['prev_lat'])
            dlambda = np.radians(self.df['lon'] - self.df['prev_lon'])
            a = np.sin(dphi/2)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda/2)**2
            c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
            self.df['dist_m'] = R * c
            y = np.sin(dlambda) * np.cos(phi2)
            x = np.cos(phi1) * np.sin(phi2) - np.sin(phi1) * np.cos(phi2) * np.cos(dlambda)
            self.df['bearing_deg'] = np.degrees(np.arctan2(y, x))

        self.df['bearing_deg'] = (self.df['bearing_deg'] + 360) % 360
        self.df['bearing_deg'] = self.df['bearing_deg'].fillna(0)
        self.df['dist_m'] = self.df['dist_m'].fillna(0)

        time_col = None
        for col in ['time', 'seconds', 'duration', 'gps_time']:
            if col in self.df.columns:
                time_col = col
                break

        if time_col:
            self.df['time_diff'] = self.df[time_col].diff()
            self.df['speed_mps'] = np.where(self.df['time_diff'] > 0,
                                            self.df['dist_m'] / self.df['time_diff'],
                                            0)
            self.df['accel_mps2'] = self.df['speed_mps'].diff().abs() / self.df['time_diff']
            self.df['accel_mps2'] = self.df['accel_mps2'].fillna(0)
        else:
            self.df['speed_mps'] = self.df['dist_m'] # fallback
            self.df['accel_mps2'] = 0.0

        self.df.drop(columns=['prev_lat', 'prev_lon'], inplace=True)

    def identify_transects_v7(self, turn_threshold=45, min_len_m=30):
        # smooth bearing (vector average)
        rads = np.radians(self.df['bearing_deg'])
        sin_avg = rads.apply(np.sin).rolling(window=3, center=True).mean()
        cos_avg = rads.apply(np.cos).rolling(window=3, center=True).mean()
        self.df['bearing_smooth'] = np.degrees(np.arctan2(sin_avg, cos_avg))
        self.df['bearing_smooth'] = (self.df['bearing_smooth'] + 360) % 360

        self.df['bearing_diff'] = self.df['bearing_smooth'].diff().abs()
        self.df['bearing_diff'] = np.where(self.df['bearing_diff'] > 180,
                                           360 - self.df['bearing_diff'],
                                           self.df['bearing_diff'])

        is_break = (self.df['bearing_diff'] > turn_threshold) | (self.df['dist_m'] > 20)
        self.df['segment_id'] = is_break.cumsum()

        stats = self.df.groupby('segment_id').agg(
            count=('lat', 'count'),
            dist_sum=('dist_m', 'sum')
        )

        valid_ids = stats[stats['dist_sum'] >= min_len_m].index
        short_ids = stats[(stats['dist_sum'] < min_len_m) & (stats['dist_sum'] > 5)].index

        transect_map = {}
        current_id = 1

        for seg_id in stats.index:
            if seg_id in valid_ids:
                transect_map[seg_id] = current_id
                current_id += 1
            elif seg_id in short_ids:
                transect_map[seg_id] = -2
            else:
                transect_map[seg_id] = -1 # Turn (very short/sharp)

        self.df['raw_tid'] = self.df['segment_id'].map(transect_map)
        self.df['transect_id'] = np.where(self.df['raw_tid'] > 0, self.df['raw_tid'], -1)
        self.df['is_turn'] = self.df['raw_tid'] == -1
        self.df['is_short'] = self.df['raw_tid'] == -2
        self.df['is_maneuver'] = (self.df['is_turn']) | (self.df['is_short'])

        self.df.drop(columns=['bearing_diff', 'segment_id', 'bearing_smooth', 'raw_tid'], inplace=True, errors='ignore')


    def identify_smart_turns(self, turn_threshold=45):
        """
        Smart Turn Detection using Pass Correlation.
        1. Identifies candidate passes using V7 logic.
        2. Correlates End of Pass(i) with Start of Pass(i+1).
        3. Analyzes geometry of the gap (Turn vs Continuation).
        4. Refines turn boundaries to include curvature entry/exit.
        """
        # Step 1: Initial rough segmentation
        self.identify_transects_v7(turn_threshold=turn_threshold)

        # Ensure we have time/index order
        if 'time' in self.df.columns:
            self.df.sort_values('time', inplace=True)

        # Get valid transects
        valid_tids = self.df[self.df['transect_id'] > 0]['transect_id'].unique()
        # sort just in case
        valid_tids = sorted(valid_tids)

        if len(valid_tids) < 2:
            return

        # Helper to get bearing diff
        def get_b_diff(b1, b2):
            d = abs(b1 - b2)
            if d > 180: d = 360 - d
            return d

        for i in range(len(valid_tids) - 1):
            t_curr = valid_tids[i]
            t_next = valid_tids[i+1]

            # Indices
            curr_indices = self.df[self.df['transect_id'] == t_curr].index
            next_indices = self.df[self.df['transect_id'] == t_next].index

            # Safety check: skip if either transect is empty
            if len(curr_indices) == 0 or len(next_indices) == 0:
                continue

            last_idx_curr = curr_indices[-1]
            first_idx_next = next_indices[0]

            # Ordering check
            if first_idx_next < last_idx_curr:
                continue

            # Gap processing
            # Rows strictly between them (by index, assuming sorted df)
            # Assuming integer index is contiguous or at least monotonic if sorted
            # But index might not be monotonic if we sorted by time and reset index?
            # Safe to use mask on boolean
            # We need to rely on the current sort order.

            # --- Correlation ---
            # Get bearing of adjacent pass ends
            N = min(5, len(curr_indices))
            b_out = self.df.loc[curr_indices[-N:], 'bearing_deg'].median()

            M = min(5, len(next_indices))
            b_in = self.df.loc[next_indices[:M], 'bearing_deg'].median()

            delta_angle = get_b_diff(b_out, b_in)

            is_u_turn = delta_angle > 135  # Approx 180

            if is_u_turn:
                # 1. Mark gap as Turn
                # Find indices between
                # If dataframe is sorted, we can use slice if we know position
                # Or just search
                # Optimization: assume indices are somewhat close
                # But safer to filter

                # Identify gap rows using index range if index is sorted?
                # Let's assume index is not necessarily sorted time-wise unless we reset it.
                # But we just sorted by time.
                # Let's relying on 'transect_id' being -1 or -2 in the gap.

                # Simple approach: Find points between the two time stamps?
                t_end_time = self.df.at[last_idx_curr, 'time'] if 'time' in self.df.columns else None
                t_start_time = self.df.at[first_idx_next, 'time'] if 'time' in self.df.columns else None

                if t_end_time is not None and t_start_time is not None:
                     gap_mask = (self.df['time'] > t_end_time) & (self.df['time'] < t_start_time)
                     self.df.loc[gap_mask, 'transect_id'] = -1
                     self.df.loc[gap_mask, 'is_turn'] = True
                     self.df.loc[gap_mask, 'is_short'] = False

                # 2. Refinement: "Eat" into the transects if they curve
                # Check End of Curr
                acc_turn = []
                for idx in reversed(curr_indices):
                    b_p = self.df.at[idx, 'bearing_deg']
                    if get_b_diff(b_p, b_out) > 10:
                        acc_turn.append(idx)
                    else:
                        break

                if acc_turn:
                    self.df.loc[acc_turn, 'transect_id'] = -1
                    self.df.loc[acc_turn, 'is_turn'] = True

                # Check Start of Next
                acc_turn = []
                for idx in next_indices:
                    b_p = self.df.at[idx, 'bearing_deg']
                    if get_b_diff(b_p, b_in) > 10:
                        acc_turn.append(idx)
                    else:
                        break

                if acc_turn:
                    self.df.loc[acc_turn, 'transect_id'] = -1
                    self.df.loc[acc_turn, 'is_turn'] = True

    def infer_swath_width(self):
        """
        Infers swath width using Total Area / Total Distance heuristic.
        Area uses CRS-aware projection (crs_utils); no preset scale.
        Fallback to 10m if fails.
        """
        for col in ['swath', 'swath_width', 'width', 'swath_ft', 'swath_m']:
            if col in self.df.columns:
                val = self.df[col].median()
                if val > 1:
                    return val

        if CRS_UTILS_AVAILABLE:
            sample = self.df if len(self.df) < 10000 else self.df.sample(10000)
            lons = sample['lon'].values
            lats = sample['lat'].values
            area_m2 = area_m2_from_points(lons, lats, self._crs)
            total_dist = self.df['dist_m'].sum()
            if total_dist > 0 and area_m2 > 0:
                calc_swath = area_m2 / total_dist
                if 3.0 < calc_swath < 30.0:
                    return calc_swath

        return 10.0

    def calculate_swath_width(self):
        return self.infer_swath_width()

    def calculate_auto_speed_limits(self):
        """
        Infers min/max speed limits from distribution.
        """
        speeds = self.df['speed_mps']
        valid_speeds = speeds[(speeds > 0.1) & (speeds < 40)] # Reasonable harvest range filter

        if valid_speeds.empty:
            return 0.5, 20.0 # Defaults

        # P05 and P99
        min_s = valid_speeds.quantile(0.05)
        max_s = valid_speeds.quantile(0.99)

        # Buffer
        min_s = max(0.1, min_s * 0.8) # 20% buffer on low end
        max_s = min(30.0, max_s * 1.2) # 20% buffer on high end

        return float(min_s), float(max_s)

    def detect_overlaps_grid_v7(self, cell_size_m=1.0):
        swath = self.calculate_swath_width()
        if self._crs and CRS_UTILS_AVAILABLE:
            x_m, y_m = points_to_meters(self.df['lon'].values, self.df['lat'].values, self._crs)
            # Center grid at origin for integer cell indices
            x_m = x_m - x_m.min()
            y_m = y_m - y_m.min()
            x_m = np.asarray(x_m, dtype=float)
            y_m = np.asarray(y_m, dtype=float)
        else:
            base_lat = self.df['lat'].min()
            base_lon = self.df['lon'].min()
            R = 6371000.0 * np.pi / 180.0
            scale_x = R * np.cos(np.radians(base_lat))
            scale_y = R
            x_m = (self.df['lon'] - base_lon).values * scale_x
            y_m = (self.df['lat'] - base_lat).values * scale_y

        grid_occupied = {}
        if 'time' in self.df.columns:
            self.df.sort_values('time', inplace=True)

        overlaps = []
        g_x = (x_m / cell_size_m).astype(int)
        g_y = (y_m / cell_size_m).astype(int)
        bearings = np.radians(self.df['bearing_deg'])

        gx_arr = np.asarray(g_x).ravel()
        gy_arr = np.asarray(g_y).ravel()
        b_arr = np.asarray(bearings).ravel()

        swath_cells = int(swath / cell_size_m)
        if swath_cells < 1: swath_cells = 1

        tids = self.df['transect_id'].values

        n_pts = len(self.df)
        report_every = max(1, n_pts // 20)  # ~5% steps
        for i in range(n_pts):
            if self._feedback and i > 0 and i % report_every == 0:
                try:
                    self._feedback.setProgress(int(25 + 25 * i / n_pts))
                    if self._feedback.isCanceled():
                        self._canceled = True
                        overlaps.extend([False] * (n_pts - i))
                        break
                except Exception:
                    pass
            bx = gx_arr[i]
            by = gy_arr[i]
            theta = b_arr[i] + (np.pi/2)
            dx = np.cos(theta) * (1.0 / cell_size_m)
            dy = np.sin(theta) * (1.0 / cell_size_m)

            # Check a short band of cells perpendicular to travel (center + +-0.25, +-0.5 swath)
            points_to_check = [(bx, by)]
            for frac in (0.25, 0.5):
                points_to_check.append((int(bx + (dx * swath_cells * frac)), int(by + (dy * swath_cells * frac))))
                points_to_check.append((int(bx - (dx * swath_cells * frac)), int(by - (dy * swath_cells * frac))))

            tid = tids[i]
            if tid < 0:
                tid = -(i + 1000000)
            is_ov = False
            for pt in points_to_check:
                if pt in grid_occupied:
                    last_tid = grid_occupied[pt]
                    if last_tid != tid:
                        is_ov = True
                else:
                    grid_occupied[pt] = tid
            overlaps.append(is_ov)

        self.df['is_overlap'] = overlaps
        return self.df['is_overlap'].sum()

    def filter_overlap_speckles(self, min_size=5):
        """
        Removes overlap flags that are not part of a sustained sequence of at least min_size.
        """
        mask = self.df['is_overlap']
        if not mask.any(): return

        # Identify changes to group consecutive values
        # True True False True True True -> groups
        change = mask.ne(mask.shift()).cumsum()

        # Count size of each group
        counts = mask.groupby(change).transform('count')

        # Filter: Keep True ONLY if count >= min_size
        # If original was False, it stays False (mask & ...)
        new_mask = mask & (counts >= min_size)

        removed_count = mask.sum() - new_mask.sum()
        if removed_count > 0:
            print(f"  [Speckle] Removed {removed_count} isolated overlap points.")

        self.df['is_overlap'] = new_mask

    def apply_filters_sequential(self, filters):
        """
        Applies filters in strict phases.
        Phase 1: Operational
        Phase 2: Global (Calculated on Op-Clean data)
        Phase 3: Spatial (Calculated on Op-Clean & Global-Clean data)
        """
        self.df['error_flag'] = 0
        var_col = filters.get('variable_col', 'yield')

        # --- PHASE 1: OPERATIONAL ---
        # 1.1 Turns & Short
        if filters.get('remove_maneuvers', True):
             self.df.loc[self.df['is_turn'], 'error_flag'] |= self.error_masks['Op: Turn']
             self.df.loc[self.df['is_short'], 'error_flag'] |= self.error_masks['Op: Short Segment']

        # 1.2 Overlaps
        self.detect_overlaps_grid_v7(cell_size_m=1.0)
        self.filter_overlap_speckles(min_size=5)
        self.df.loc[self.df['is_overlap'], 'error_flag'] |= self.error_masks['Op: Overlap']

        # 1.3 Speed & Accel
        if 'min_speed' in filters:
            self.df.loc[self.df['speed_mps'] < filters['min_speed'], 'error_flag'] |= self.error_masks['Op: Min Speed']
        if 'max_speed' in filters:
            self.df.loc[self.df['speed_mps'] > filters['max_speed'], 'error_flag'] |= self.error_masks['Op: Max Speed']
        if 'accel_mps2' in self.df.columns:
            self.df.loc[self.df['accel_mps2'] > 1.5, 'error_flag'] |= self.error_masks['Op: High Acceleration']

        # 1.4 Distance Jumps
        median_dist = self.df['dist_m'].median()
        if median_dist > 0:
            self.df.loc[self.df['dist_m'] > (5 * median_dist), 'error_flag'] |= self.error_masks['Op: Distance Jump']

        # 1.5 Smooth Speed (Inspired by YE Smooth Velocity)
        if filters.get('speed_smooth_factor', 0) > 0:
            self.df['speed_diff_ratio'] = self.df['speed_mps'].diff().abs() / self.df['speed_mps']
            self.df.loc[self.df['speed_diff_ratio'] > filters['speed_smooth_factor'], 'error_flag'] |= self.error_masks['Op: Variable Speed Change']

        # 1.6 Delays
        # Apply to points NOT YET FLAGGED (or apply to all and merge flags? Sequential implies clean baseline for detection?)
        # For Delay, we look at Transect structure. Maneuvers are excluded.
        # But if speed is low, it might affect delay calc.
        # Let's clean the DataFrame temporarily to find delays?
        # Creating a 'clean_for_delay' mask
        # Actually delay detection logic usually ignores maneuvers.
        auto_start = filters.get('auto_start_delay', False)
        auto_end = filters.get('auto_end_delay', False)
        if (auto_start or auto_end) and var_col in self.df.columns:
            self.apply_per_transect_delay(var_col, enable_start=auto_start, enable_end=auto_end)

        # --- PHASE 2: GLOBAL ---
        # Only evaluate on Operational-Clean points
        op_mask = (self.df['error_flag'] == 0)

        if var_col in self.df.columns:
            val_series = self.df.loc[op_mask, var_col]

            # Limits
            min_y = filters.get('min_yield', None)
            max_y = filters.get('max_yield', None)

            # Auto Limits (AOI)
            if filters.get('use_aoi', False) and not val_series.empty:
                low_a, high_a = self.calculate_aoi_ranges(val_series)
                # User-adjustable margin (default 0.3 = 30% wider bounds)
                margin = float(filters.get('global_outlier_margin', 0.3))
                margin = max(0.0, min(1.0, margin))
                low_a = low_a - low_a * margin
                high_a = high_a + high_a * margin

                # If explicit min not set, use AOI
                if min_y is None: min_y = low_a
                else: min_y = max(min_y, low_a) # Stricter? or Union? Sequential = user overrides?
                # Usually User Limits > AOI auto.

                if max_y is None: max_y = high_a
                else: max_y = min(max_y, high_a)

            # Apply Limit Checks only to CLEAN points?
            # User says: "If first wave takes out a point, it is no longer considered"
            # So we only flag points that are currently Clean.
            if min_y is not None:
                self.df.loc[op_mask & (self.df[var_col] < min_y), 'error_flag'] |= self.error_masks['Gl: Value Low']
            if max_y is not None:
                self.df.loc[op_mask & (self.df[var_col] > max_y), 'error_flag'] |= self.error_masks['Gl: Value High']


        # --- PHASE 3: SPATIAL ---
        # Only use strictly CLEAN points (Op + Global Clean) for stats
        # Flagging can apply to currently Clean points.

        if filters.get('use_local_std', True) and var_col in self.df.columns:
            sigma = filters.get('local_std_sigma', 3.0)
            cell_mult = filters.get('grid_cell_multiplier', 3.0)
            self.apply_local_std_filter(var_col, sigma, cell_mult)

        if filters.get('use_cross_track', False) or filters.get('use_cross_track_correct', False):
            self.apply_cross_track_filter(var_col, filters)

        # Labels
        self.df['error_reason'] = self.df['error_flag'].apply(lambda x: "; ".join([k for k,v in self.error_masks.items() if x & v]) if x > 0 else "Clean")

        op_mask_val = 1 | 2 | 4 | 8 | 16 | 32 | 64 | 1024 | 2048 | 8192
        gl_mask_val = 128 | 256
        lo_mask_val = 512 | 4096

        def get_category(flag):
            if flag == 0: return "Clean"
            if flag & op_mask_val: return "Operational Error"
            if flag & gl_mask_val: return "Global Outlier"
            if flag & lo_mask_val: return "Local Outlier"
            return "Unknown"

        self.df['filtering_category'] = self.df['error_flag'].apply(get_category)

        # Optional: simple annotation (Clean, Local, Global, Operational Error only)
        if filters.get('simple_annotation', False):
            self._apply_simple_annotation()

        # Output column with same values as filtering_category (for display as "Error Type")
        self.df['Error Type'] = self.df['filtering_category']

    def _apply_simple_annotation(self):
        """
        Simple 4-way labels, sequential: Phase 1 Operational -> Phase 2 Global -> Phase 3 Local -> Clean.
        Phase 1: any point with transect_id < 0 (turn/short) OR operational bits in error_flag -> "Operational Error".
        """
        gl_mask = 128 | 256
        lo_mask = 512 | 4096   # Spatial Outlier, Cross-Track (Overlap 64 -> Operational Error)
        op_flag_mask = 1 | 2 | 16 | 32 | 64 | 1024 | 2048   # + Overlap (64)

        def safe_flag(flag):
            try:
                if flag is None or (isinstance(flag, float) and np.isnan(flag)):
                    return 0
                return int(flag)
            except (TypeError, ValueError):
                return 0

        n = len(self.df)
        simple = np.full(n, "Clean", dtype=object)
        flags = self.df["error_flag"].map(safe_flag).values

        # Phase 1: Operational Error — use every available signal so turns are never missed
        maneuver = np.zeros(n, dtype=bool)
        if "transect_id" in self.df.columns:
            tid = np.asarray(self.df["transect_id"].values, dtype=float)
            maneuver = (tid < 0) & ~np.isnan(tid)
        if "is_maneuver" in self.df.columns:
            maneuver = maneuver | np.asarray(self.df["is_maneuver"].fillna(False).values, dtype=bool)
        if "is_turn" in self.df.columns:
            maneuver = maneuver | np.asarray(self.df["is_turn"].fillna(False).values, dtype=bool)
        if "is_short" in self.df.columns:
            maneuver = maneuver | np.asarray(self.df["is_short"].fillna(False).values, dtype=bool)
        op_from_flags = (flags & op_flag_mask) != 0
        op_hit = maneuver | op_from_flags
        simple[op_hit] = "Operational Error"

        # Phase 2: Global (only where still Clean)
        still_clean = simple == "Clean"
        simple[still_clean & ((flags & gl_mask) != 0)] = "Global"

        # Phase 3: Local only for points IN a valid harvest pass (transect_id > 0).
        # Points not in a valid pass (transect_id <= 0) with local flags -> Operational Error, never Local.
        still_clean = simple == "Clean"
        has_local_flag = (flags & lo_mask) != 0
        in_valid_pass = np.ones(n, dtype=bool)
        if "transect_id" in self.df.columns:
            tid = np.asarray(self.df["transect_id"].values, dtype=float)
            in_valid_pass = (tid > 0) & ~np.isnan(tid)
        assign_local = still_clean & has_local_flag & in_valid_pass
        assign_op_instead = still_clean & has_local_flag & ~in_valid_pass  # turn zone but had local flag
        simple[assign_local] = "Local"
        simple[assign_op_instead] = "Operational Error"

        # Safety: any point that is turn/short by geometry must be Operational Error, not Local
        if maneuver.any():
            wrong_local = (simple == "Local") & maneuver
            simple[wrong_local] = "Operational Error"

        self.df["error_reason"] = simple
        self.df["filtering_category"] = simple

    def medcouple(self, x):
        x = np.array(x, dtype=np.float64)
        x = x[~np.isnan(x)]
        n = len(x)
        if n < 3: return 0.0



        x_sorted = np.sort(x)
        median_x = np.median(x_sorted)
        z = x_sorted - median_x
        lower = z[z <= 0]
        upper = z[z >= 0]

        if len(lower) == 0 or len(upper) == 0: return 0.0

        if n > 2000:
            np.random.seed(42)
            indices = np.random.choice(n, size=2000, replace=False)
            x_sub = x[indices]
            return self.medcouple(x_sub)

        lower_g, upper_g = np.meshgrid(lower, upper)
        diff = upper_g - lower_g
        sum_vals = upper_g + lower_g

        with np.errstate(divide='ignore', invalid='ignore'):
            kernel = sum_vals / diff
        kernel[diff == 0] = np.sign(sum_vals[diff == 0])
        return np.median(kernel)

    def calculate_aoi_ranges(self, data, k=1.5):
        series = pd.Series(data).dropna()
        if len(series) < 5: return -np.inf, np.inf
        Q1 = series.quantile(0.25)
        Q3 = series.quantile(0.75)
        IQR = Q3 - Q1
        MC = self.medcouple(series.values)
        if MC >= 0:
            lower = Q1 - k * np.exp(-4 * MC) * IQR
            upper = Q3 + k * np.exp(3 * MC) * IQR
        else:
            lower = Q1 - k * np.exp(-3 * MC) * IQR
            upper = Q3 + k * np.exp(4 * MC) * IQR
        return float(lower), float(upper)

    def apply_local_std_filter(self, var_col, sigma, cell_multiplier):
        swath = self.calculate_swath_width()
        cell_size_m = swath * cell_multiplier
        if cell_size_m <= 0:
            cell_size_m = 0.1
        if self._crs and CRS_UTILS_AVAILABLE:
            x_m, y_m = points_to_meters(self.df['lon'].values, self.df['lat'].values, self._crs)
            self.df['lx'] = (np.asarray(x_m, dtype=float) / cell_size_m).astype(int)
            self.df['ly'] = (np.asarray(y_m, dtype=float) / cell_size_m).astype(int)
        else:
            # Fallback: degree-based cell size (approx at mid-lat)
            R = 6371000.0 * np.pi / 180.0
            center_lat = self.df['lat'].mean()
            deg_per_m = 1.0 / (R * np.cos(np.radians(center_lat)))
            cell_deg = cell_size_m * deg_per_m
            if cell_deg <= 0:
                cell_deg = 1e-6
            self.df['lx'] = (self.df['lon'] / cell_deg).astype(int)
            self.df['ly'] = (self.df['lat'] / cell_deg).astype(int)

        # SEQUENTIAL: Use only currently clean data for stats
        mask_clean = (self.df['error_flag'] == 0) & (self.df[var_col].notna())

        stats = self.df[mask_clean].groupby(['lx', 'ly'])[var_col].agg(['mean', 'std'])
        self.df = self.df.merge(stats, on=['lx', 'ly'], how='left')

        self.df['std'] = self.df['std'].fillna(0)
        abnormal = (self.df[var_col] - self.df['mean']).abs() > (sigma * self.df['std'])
        abnormal = abnormal & (self.df['std'] > 0)

        # Only flag if currently clean (Sequential)
        self.df.loc[abnormal & mask_clean, 'error_flag'] |= self.error_masks['Lo: Spatial Outlier']
        self.df.drop(columns=['lx', 'ly', 'mean', 'std'], inplace=True, errors='ignore')

    def apply_cross_track_filter(self, var_col, filters=None):
        if 'transect_id' not in self.df.columns:
            return
        filters = filters or {}
        clean_df = self.df[self.df['error_flag'] == 0]
        t_medians = clean_df[clean_df['transect_id'] > 0].groupby('transect_id')[var_col].median()
        if t_medians.empty:
            return
        t_roll = t_medians.rolling(window=5, center=True).median()
        ratio = t_medians / t_roll.replace(0, np.nan)
        bad_tids = ratio[(ratio > 1.5) | (ratio < 0.66)].index

        # Flag stripy transects (optional)
        if filters.get('use_cross_track', False):
            self.df.loc[self.df['transect_id'].isin(bad_tids) & (self.df['error_flag'] == 0), 'error_flag'] |= self.error_masks['Lo: Cross-Track Outlier']

        # Correct stripy transects (optional): scale each transect to neighborhood median
        if filters.get('use_cross_track_correct', False):
            factors = (t_roll / t_medians).replace([np.inf, -np.inf, np.nan], 1.0)
            factors = factors.clip(0.1, 10.0)
            self.df[var_col + '_corrected'] = (
                self.df[var_col].astype(float) * self.df['transect_id'].map(factors).fillna(1.0)
            )

    def apply_per_transect_delay(self, var_col, enable_start=True, enable_end=True):
         if 'transect_id' not in self.df.columns: return

         valid_tids = self.df[self.df['transect_id'] > 0]['transect_id'].unique()
         start_indices = []
         end_indices = []

         for tid in valid_tids:
            rows = self.df[self.df['transect_id'] == tid]
            if len(rows) < 5: continue

            # Use median of the transect (robust to outliers, but ideally should use Clean data?)
            # If we haven't filtered Global yet, we might have outliers.
            # But delay depends on ramp. Robust median is usually fine.
            median_val = rows[var_col].median()
            if median_val <= 0: continue

            thresh = 0.70 * median_val
            vals = rows[var_col].values
            indices = rows.index.values

            if enable_start:
                limit = min(len(vals), 20)
                cut = 0
                for i in range(limit):
                    if vals[i] >= thresh:
                        cut = i
                        break
                if cut > 0:
                    start_indices.extend(indices[:cut])

            if enable_end:
                limit = min(len(vals), 20)
                cut = 0
                N = len(vals)
                for i in range(limit):
                    if vals[N - 1 - i] >= thresh:
                        cut = i
                        break
                if cut > 0:
                    end_indices.extend(indices[N - cut:])

         if start_indices:
            self.df.loc[start_indices, 'error_flag'] |= self.error_masks['Op: Start Delay']
         if end_indices:
            self.df.loc[end_indices, 'error_flag'] |= self.error_masks['Op: End Delay']

    def infer_delay_parameters(self, var_col):
        """
        Analyzes the first 20 points of passes to see if there is a 'Ramp Up'.
        Returns (enable_start_delay, enable_end_delay)
        """
        if 'transect_id' not in self.df.columns:
             return False, False

        valid_tids = self.df[self.df['transect_id'] > 0]['transect_id'].unique()
        if len(valid_tids) < 5: return False, False

        # Sample random transects
        sample_size = min(len(valid_tids), 30)
        sample_tids = np.random.choice(valid_tids, size=sample_size, replace=False)

        ramp_count = 0
        total_checked = 0

        for tid in sample_tids:
            rows = self.df[self.df['transect_id'] == tid]
            if len(rows) < 20: continue

            # Check Start Ramp
            # Logic: Correlation between Index (0..10) and Yield should be positive and high
            vals = rows[var_col].iloc[:10].values
            if len(vals) < 5: continue

            # Simple check: Is mean of first 5 < mean of next 5?
            m1 = np.mean(vals[:5])
            m2 = np.mean(vals[5:])

            # Significant ramp? > 20% growth
            if m2 > (m1 * 1.2):
                ramp_count += 1
            total_checked += 1

        if total_checked == 0: return False, False

        # If > 50% have ramp, enable it
        enable_start = (ramp_count / total_checked) > 0.5

        return enable_start, False # Only auto-detect start for now

    def get_cleaned_df(self):
        return self.df

    def get_stats(self, var_col):
        return {
            'min_var': self.df[var_col].min(), 'max_var': self.df[var_col].max(),
            'min_speed': self.df['speed_mps'].min(), 'max_speed': self.df['speed_mps'].max()
        }


# Alias for backward compatibility: test_smart_turn.py, test_automation.py,
# and verify_import.py import YieldCleaner from this module.
YieldCleaner = DataCleaner


# Processing Algorithm Wrapper
class SpatialDataFilteringAlgorithm(QgsProcessingAlgorithm):
    INPUT_LAYER = 'INPUT_LAYER'
    VARIABLE_FIELD = 'VARIABLE_FIELD'
    AUTO_START_DELAY = 'AUTO_START_DELAY'
    AUTO_END_DELAY = 'AUTO_END_DELAY'
    MIN_SPEED = 'MIN_SPEED'
    MAX_SPEED = 'MAX_SPEED'
    # New Global Params
    MIN_VALUE = 'MIN_VALUE'
    MAX_VALUE = 'MAX_VALUE'
    USE_AOI = 'USE_AOI'
    GLOBAL_OUTLIER_MARGIN = 'GLOBAL_OUTLIER_MARGIN'

    TURN_THRESHOLD_DEG = 'TURN_THRESHOLD_DEG'
    AUTO_DETECT_TURN = 'AUTO_DETECT_TURN'

    LOCAL_STD_SIGMA = 'LOCAL_STD_SIGMA'
    GRID_CELL_MULTIPLIER = 'GRID_CELL_MULTIPLIER'
    USE_CROSS_TRACK = 'USE_CROSS_TRACK'
    USE_CROSS_TRACK_CORRECT = 'USE_CROSS_TRACK_CORRECT'
    SIMPLE_ANNOTATION = 'SIMPLE_ANNOTATION'
    OUTPUT = 'OUTPUT'

    def createInstance(self): return SpatialDataFilteringAlgorithm()
    def name(self): return 'spatialdatafiltering'
    def displayName(self): return 'Spatial Data Filtering (V8 Sequential)'
    def group(self): return 'RETA'
    def groupId(self): return 'reta'

    AUTO_SPEED = 'AUTO_SPEED'

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT_LAYER, 'Input Layer', [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterField(self.VARIABLE_FIELD, 'Variable to Clean', parentLayerParameterName=self.INPUT_LAYER, type=QgsProcessingParameterField.Numeric))

        # --- AUTOMATION SECTION ---
        self.addParameter(QgsProcessingParameterBoolean(self.AUTO_SPEED, 'Auto-Detect Speeds', defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(self.AUTO_START_DELAY, 'Auto-Detect Start Delays', defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(self.AUTO_END_DELAY, 'Auto-Detect End Delays', defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(self.AUTO_DETECT_TURN, 'Auto-Detect Turn Threshold', defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(self.USE_AOI, 'Auto-Yield Limits (AOI)', defaultValue=True))
        p_margin = QgsProcessingParameterNumber(
            self.GLOBAL_OUTLIER_MARGIN,
            'Global outlier margin (0-1; fraction to widen AOI bounds; 0.3 = 30%)',
            type=QgsProcessingParameterNumber.Double,
            defaultValue=0.3,
            optional=True
        )
        p_margin.setFlags(p_margin.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(p_margin)

        # --- MANUAL OVERRIDES (ADVANCED) ---
        # Note: We use FlagAdvanced (1) to hide them by default

        # 1. Speed
        p_min_s = QgsProcessingParameterNumber(self.MIN_SPEED, 'Min Speed (m/s) [Ignored if Auto-Speed True]', type=QgsProcessingParameterNumber.Double, defaultValue=0.2)
        p_min_s.setFlags(p_min_s.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(p_min_s)

        p_max_s = QgsProcessingParameterNumber(self.MAX_SPEED, 'Max Speed (m/s) [Ignored if Auto-Speed True]', type=QgsProcessingParameterNumber.Double, defaultValue=20.0)
        p_max_s.setFlags(p_max_s.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(p_max_s)

        # 2. Yield Limits
        p_min_y = QgsProcessingParameterNumber(self.MIN_VALUE, 'Min Value [Ignored if Auto-Yield True]', type=QgsProcessingParameterNumber.Double, optional=True, defaultValue=0.0)
        p_min_y.setFlags(p_min_y.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(p_min_y)

        p_max_y = QgsProcessingParameterNumber(self.MAX_VALUE, 'Max Value [Ignored if Auto-Yield True]', type=QgsProcessingParameterNumber.Double, optional=True, defaultValue=2000.0)
        p_max_y.setFlags(p_max_y.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(p_max_y)

        # 3. Turn Threshold
        p_turn = QgsProcessingParameterNumber(self.TURN_THRESHOLD_DEG, 'Turn Threshold (Degrees)', type=QgsProcessingParameterNumber.Double, defaultValue=15.0)
        p_turn.setFlags(p_turn.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(p_turn)

        # Other
        self.addParameter(QgsProcessingParameterNumber(self.LOCAL_STD_SIGMA, 'Local Outlier Sigma', type=QgsProcessingParameterNumber.Double, defaultValue=3.0))
        self.addParameter(QgsProcessingParameterNumber(self.GRID_CELL_MULTIPLIER, 'Grid Cell Multiplier (x Swath)', type=QgsProcessingParameterNumber.Double, defaultValue=3.0))
        self.addParameter(QgsProcessingParameterBoolean(self.USE_CROSS_TRACK, 'Flag stripy transects (Cross-Track)', defaultValue=True))
        self.addParameter(QgsProcessingParameterBoolean(self.USE_CROSS_TRACK_CORRECT, 'Correct stripy transects (write corrected variable)', defaultValue=False))
        self.addParameter(QgsProcessingParameterBoolean(self.SIMPLE_ANNOTATION, 'Simple error categories (Clean, Local, Global, Operational Error only)', defaultValue=False))
        self.addParameter(QgsProcessingParameterFeatureSink(self.OUTPUT, 'Filtered Output'))

    def processAlgorithm(self, parameters, context, feedback):
        try:
            feedback.setProgress(0)
        except Exception:
            pass
        source = self.parameterAsSource(parameters, self.INPUT_LAYER, context)
        var_field = self.parameterAsString(parameters, self.VARIABLE_FIELD, context)
        turn_thresh = self.parameterAsDouble(parameters, self.TURN_THRESHOLD_DEG, context)

        features = []
        fields = [f.name() for f in source.fields()]
        for feat in source.getFeatures():
            dct = dict(zip(fields, feat.attributes()))
            geom = feat.geometry()
            if not geom.isEmpty():
                dct['lat'] = geom.asPoint().y()
                dct['lon'] = geom.asPoint().x()
                features.append(dct)

        if not features:
            return {self.OUTPUT: None}
        df = pd.DataFrame(features)
        try:
            feedback.setProgress(5)
        except Exception:
            pass

        # CRS: use layer CRS for proper distance/area (no preset scale)
        crs = None
        try:
            if source.sourceCrs().isValid():
                crs = source.sourceCrs()
                feedback.pushInfo(f"[CRS] Using layer CRS: {crs.authid()}")
        except Exception:
            pass
        cleaner = DataCleaner(df, crs=crs, feedback=feedback)
        cleaner.calculate_geometry()
        try:
            feedback.setProgress(10)
            if feedback.isCanceled():
                return {self.OUTPUT: None}
        except Exception:
            pass

        # --- AUTO INFERENCE (always run when auto is on; user can override via manual params) ---

        # 1. Turn threshold
        auto_turn = self.parameterAsBool(parameters, self.AUTO_DETECT_TURN, context)
        if auto_turn:
            feedback.pushInfo("[Auto] Detecting Smart Turn Threshold...")
            cleaner.identify_smart_turns(turn_threshold=turn_thresh)
            feedback.pushInfo(f"[Auto] Turn threshold in use: {turn_thresh} (uncheck Auto to set manually)")
        else:
            cleaner.identify_transects_v7(turn_threshold=turn_thresh)
            feedback.pushInfo(f"[Manual] Turn threshold: {turn_thresh}")
        try:
            feedback.setProgress(15)
            if feedback.isCanceled():
                return {self.OUTPUT: None}
        except Exception:
            pass

        # 2. Speed Limits
        auto_speed = self.parameterAsBool(parameters, self.AUTO_SPEED, context)

        if auto_speed:
            auto_min_s, auto_max_s = cleaner.calculate_auto_speed_limits()
            feedback.pushInfo(f"[Auto] Inferred speed limits: {auto_min_s:.2f}-{auto_max_s:.2f} m/s (uncheck Auto to type your own)")
            min_s_final = auto_min_s
            max_s_final = auto_max_s
        else:
            min_s_final = self.parameterAsDouble(parameters, self.MIN_SPEED, context)
            max_s_final = self.parameterAsDouble(parameters, self.MAX_SPEED, context)
            feedback.pushInfo(f"[Manual] Speed limits: {min_s_final:.2f}-{max_s_final:.2f} m/s")
        try:
            feedback.setProgress(20)
            if feedback.isCanceled():
                return {self.OUTPUT: None}
        except Exception:
            pass

        # 3. Delay Detection
        p_start_delay = self.parameterAsBool(parameters, self.AUTO_START_DELAY, context)

        inferred_start, _ = cleaner.infer_delay_parameters(var_field)

        # Logic: If Auto-Start is Checked, we use Inference.
        # But wait, if user UNCHECKS it, they might want to force it OFF?
        # Or force it ON?
        # Standard QGIS Boolean is just a Flag.
        # Interpretation: "Enable Auto-Detection". If True -> Run Inference. If Inference=True -> Apply.
        # If False -> Do not run inference (Assume Off).
        # This matches "Auto" behavior.

        use_start_delay = False
        if p_start_delay:
            if inferred_start:
                feedback.pushInfo("[Auto] Detected start delay signature; enabling correction.")
                use_start_delay = True
            else:
                feedback.pushInfo("[Auto] No start delay signature; leaving off.")
        else:
            feedback.pushInfo("[Manual] Start delay: off (check Auto to detect)")

        feedback.pushInfo(f"Stats: {cleaner.get_stats(var_field)}")

        # 4. Yield limits (AOI) and global outlier margin
        use_aoi = self.parameterAsBool(parameters, self.USE_AOI, context)
        min_y_param = self.parameterAsDouble(parameters, self.MIN_VALUE, context)
        max_y_param = self.parameterAsDouble(parameters, self.MAX_VALUE, context)
        try:
            margin_param = self.parameterAsDouble(parameters, self.GLOBAL_OUTLIER_MARGIN, context)
        except Exception:
            margin_param = 0.3
        if use_aoi and len(cleaner.df) > 0 and var_field in cleaner.df.columns:
            low_a, high_a = cleaner.calculate_aoi_ranges(cleaner.df[var_field].dropna())
            margin = max(0.0, min(1.0, margin_param))
            low_display = low_a - low_a * margin
            high_display = high_a + high_a * margin
            feedback.pushInfo(f"[Auto] AOI bounds (margin={margin:.2f}): [{low_display:.2f}, {high_display:.2f}] (uncheck Auto to set min/max manually)")
        elif not use_aoi:
            feedback.pushInfo(f"[Manual] Yield limits: min={min_y_param}, max={max_y_param}")

        filters = {
            'variable_col': var_field,
            'auto_start_delay': use_start_delay,
            'auto_end_delay': self.parameterAsBool(parameters, self.AUTO_END_DELAY, context),
            'min_speed': min_s_final,
            'max_speed': max_s_final,
            'min_yield': min_y_param,
            'max_yield': max_y_param,
            'use_aoi': use_aoi,
            'global_outlier_margin': margin_param,
            'local_std_sigma': self.parameterAsDouble(parameters, self.LOCAL_STD_SIGMA, context),
            'grid_cell_multiplier': self.parameterAsDouble(parameters, self.GRID_CELL_MULTIPLIER, context),
            'use_cross_track': self.parameterAsBool(parameters, self.USE_CROSS_TRACK, context),
            'use_cross_track_correct': self.parameterAsBool(parameters, self.USE_CROSS_TRACK_CORRECT, context),
            'simple_annotation': self.parameterAsBool(parameters, self.SIMPLE_ANNOTATION, context),
            'remove_maneuvers': True
        }
        try:
            feedback.setProgress(25)
        except Exception:
            pass
        cleaner.apply_filters_sequential(filters)
        if getattr(cleaner, '_canceled', False):
            return {self.OUTPUT: None}
        try:
            feedback.setProgress(85)
        except Exception:
            pass
        res = cleaner.get_cleaned_df()
        use_cross_track_correct = self.parameterAsBool(parameters, self.USE_CROSS_TRACK_CORRECT, context)
        corrected_col = var_field + '_corrected' if use_cross_track_correct else None
        if corrected_col and corrected_col not in res.columns:
            corrected_col = None

        out_fields = source.fields()
        new_cols = ['dist_m', 'bearing_deg', 'speed_mps', 'transect_id', 'is_turn', 'is_short', 'error_flag', 'error_reason', 'filtering_category', 'Error Type']
        if corrected_col:
            new_cols.append(corrected_col)
        for c in new_cols:
            t = QVariant.String if c in ['error_reason', 'filtering_category', 'Error Type'] else (QVariant.Bool if 'is_' in c else QVariant.Double)
            if 'flag' in c or 'id' in c: t = QVariant.Int
            out_fields.append(QgsField(c, t))

        sink, dest_id = self.parameterAsSink(parameters, self.OUTPUT, context, out_fields, source.wkbType(), source.sourceCrs())
        n_out = len(res)
        for idx, (_, row) in enumerate(res.iterrows()):
            try:
                if idx > 0 and idx % max(1, n_out // 15) == 0:
                    feedback.setProgress(85 + int(15 * idx / n_out))
                    if feedback.isCanceled():
                        return {self.OUTPUT: None}
            except Exception:
                pass
            fet = QgsFeature()
            fet.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(row['lon'], row['lat'])))
            attrs = [row[f.name()] if f.name() in row else None for f in source.fields()]
            attrs.append(float(row.get('dist_m', 0)))
            attrs.append(float(row.get('bearing_deg', 0)))
            attrs.append(float(row.get('speed_mps', 0)))
            attrs.append(int(row.get('transect_id', -1)))
            attrs.append(bool(row.get('is_turn', False)))
            attrs.append(bool(row.get('is_short', False)))
            attrs.append(int(row.get('error_flag', 0)))
            attrs.append(str(row.get('error_reason', '')))
            attrs.append(str(row.get('filtering_category', '')))
            attrs.append(str(row.get('Error Type', '')))
            if corrected_col:
                v = row.get(corrected_col)
                attrs.append(float(v) if v is not None else None)
            fet.setAttributes(attrs)
            sink.addFeature(fet, QgsFeatureSink.FastInsert)
        try:
            feedback.setProgress(100)
        except Exception:
            pass
        return {self.OUTPUT: dest_id}
