<p align="center">
  <img src="RETA_logo.png" alt="RETA Logo" width="200"/>
</p>

<h1 align="center">RETA - Robust Error-treatment Toolkit for Agriculture</h1>

<p align="center">
  A QGIS Processing plugin for sequential spatial data filtering of point-based sensor data.
</p>

---

## Overview

RETA provides a three-phase sequential filtering pipeline for cleaning spatial point data collected by mobile platforms (e.g., combine harvesters, soil sensors, electromagnetic induction instruments). It flags erroneous observations using bitmask error codes and categorizes each point as **Clean**, **Operational Error**, **Global Outlier**, or **Local Outlier**.

While applicable to any georeferenced point dataset with a numeric variable of interest, RETA was originally adapted from and builds upon the methods described in **Yield Editor 2.0**:

> Sudduth, K. A., Drummond, S. T., & Myers, D. B. (2012, July 29--August 1). *Yield Editor 2.0: Software for automated removal of yield map errors* [Paper presentation]. ASABE Annual International Meeting, Dallas, TX, United States. https://doi.org/10.13031/2013.41893

## Authors

- **Zachary Komarnisky** and **Felippe Karp**
- Contact: zachary.komarnisky@oldscollege.ca

Agentic coding AI was used to develop this tool with human guidance and oversight.

## Installation

1. Download or clone this repository.
2. Copy the folder into your QGIS plugins directory:
   - **Windows:** `%APPDATA%/QGIS/QGIS3/profiles/default/python/plugins/`
   - **Linux:** `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
   - **macOS:** `~/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins/`
3. Restart QGIS and enable the plugin in **Plugins > Manage and Install Plugins**.
4. The algorithm appears in the **Processing Toolbox** under **RETA > Spatial Data Filtering**.

### Dependencies

RETA uses libraries bundled with QGIS. No additional installation is required in most cases.

| Library | Usage | License |
|---|---|---|
| numpy | Array operations, trigonometry | BSD-3-Clause |
| pandas | DataFrame manipulation | BSD-3-Clause |
| pyproj | CRS-aware coordinate transforms | MIT |
| scipy | Convex hull for area estimation (optional) | BSD-3-Clause |

## Filtering Pipeline

RETA applies filters in three strict sequential phases. Each subsequent phase only evaluates points that passed all previous phases, preventing cascading errors.

```
Input Layer
    |
    v
[Geometry Calculation] -- distance, bearing, speed, acceleration
    |
    v
[Transect Identification] -- segment passes, detect turns
    |
    v
Phase 1: OPERATIONAL FILTERS
    |  Turn/Short segment removal
    |  Overlap detection (grid-based)
    |  Speed limits (min/max)
    |  High acceleration
    |  Distance jumps
    |  Start/End delay (sensor-specific)
    |  Smooth velocity change
    |
    v
Phase 2: GLOBAL OUTLIER FILTERS (on Phase 1-clean data only)
    |  Adjusted boxplot (AOI) using medcouple
    |  Manual min/max value thresholds
    |
    v
Phase 3: LOCAL/SPATIAL FILTERS (on Phase 1+2-clean data only)
    |  Grid-based local standard deviation
    |  Cross-track (stripy transect) detection/correction
    |
    v
Output Layer (with error_flag, error_reason, filtering_category)
```

## Error Taxonomy

Each point receives a bitmask `error_flag` that encodes exactly which filters flagged it:

| Flag | Bit | Category | Description |
|---|---|---|---|
| `Clean` | 0 | -- | No errors detected |
| `Op: Turn` | 1 | Operational | Point is in a turn/headland maneuver |
| `Op: Short Segment` | 2 | Operational | Point is in a segment shorter than minimum length |
| `Op: Distance Jump` | 4 | Operational | Distance to previous point exceeds 5x the median |
| `Op: Collocated` | 8 | Operational | Co-located with another point |
| `Op: Start Delay` | 16 | Operational | Sensor ramp-up at pass start |
| `Op: End Delay` | 32 | Operational | Sensor ramp-down at pass end |
| `Op: Overlap` | 64 | Operational | Overlapping with a previous pass |
| `Gl: Value Low` | 128 | Global | Below global lower bound |
| `Gl: Value High` | 256 | Global | Above global upper bound |
| `Lo: Spatial Outlier` | 512 | Local | Deviates from local neighborhood |
| `Op: Min Speed` | 1024 | Operational | Below minimum speed threshold |
| `Op: Max Speed` | 2048 | Operational | Above maximum speed threshold |
| `Lo: Cross-Track Outlier` | 4096 | Local | Transect median deviates from neighbors |
| `Op: High Acceleration` | 8192 | Operational | Acceleration exceeds threshold |
| `Op: Variable Speed Change` | 16384 | Operational | Rapid speed change between consecutive points |

---

## Filter Details and Equations

### Pre-processing: Geometry Calculation

Coordinates are projected to a local metric system (UTM via pyproj when available, or approximate degree-to-meter scaling as fallback). All subsequent distance and area calculations use these projected coordinates.

**Distance** between consecutive points on the projected plane:

$$d_i = \sqrt{(\Delta x_i)^2 + (\Delta y_i)^2}$$

**Bearing** (compass convention, clockwise from north):

$$\theta_i = \left(90° - \arctan2(\Delta y_i,\, \Delta x_i) \cdot \frac{180}{\pi}\right) \bmod 360$$

**Speed** and **acceleration** (when timestamps are available):

$$v_i = \frac{d_i}{\Delta t_i} \qquad a_i = \frac{|v_i - v_{i-1}|}{\Delta t_i}$$

### Pre-processing: Transect Identification

Bearing is smoothed with a centered rolling vector average (window = 3). Segments are split where the smoothed bearing change exceeds a threshold or the inter-point distance exceeds 20 m. Segments are classified as:

- **Valid pass**: cumulative length >= `min_len_m` (default 30 m)
- **Short segment**: length between 5 m and `min_len_m`
- **Turn**: length < 5 m or sharp bearing change

**Smart Turn Detection** refines boundaries by correlating the exit bearing of pass *i* with the entry bearing of pass *i+1*. If the angular difference exceeds 135 deg (approximate U-turn), the gap is marked as a turn and boundary points showing curvature (> 10 deg deviation from the pass median bearing) are absorbed into the turn zone.

---

### Phase 1: Operational Filters

#### 1.1 Turn and Short Segment Removal

Points identified as turns or short segments during transect identification are flagged. These represent headland maneuvers, direction changes, or segments too short to contain reliable data.

#### 1.2 Overlap Detection (Grid Occupancy)

A hash-map of 1 m grid cells tracks which transect first occupied each cell. For each point, a band of cells perpendicular to the travel direction (at fractions of the swath width) is checked. If any cell was previously occupied by a different transect, the point is flagged as overlap.

Isolated overlap flags (runs shorter than 5 consecutive points) are removed as speckle noise.

*Reference: Sudduth & Drummond (2007) describe overlap filtering based on swath geometry and grid-cell occupancy.*

> Sudduth, K. A., & Drummond, S. T. (2007). Yield Editor: Software for removing errors from crop yield maps. *Agronomy Journal*, *99*(6), 1471--1482. https://doi.org/10.2134/agronj2006.0326

#### 1.3 Speed Filtering

Points where platform speed falls outside acceptable bounds are flagged:

$$\text{Flag if: } v_i < v_{\min} \quad \text{or} \quad v_i > v_{\max}$$

In auto-detect mode, bounds are derived from the 5th and 99th percentiles of the speed distribution with a 20% buffer.

*This is particularly important for yield data, where speed approaching zero causes computed yield values to become artificially inflated (mass flow / speed). Adapted from the MINV and MAXV filters in Yield Editor (Sudduth & Drummond, 2007).*

#### 1.4 High Acceleration

Points with acceleration exceeding 1.5 m/s^2 are flagged, indicating abrupt changes in platform motion that may compromise sensor readings.

#### 1.5 Distance Jumps

Points where the inter-point distance exceeds 5 times the dataset median distance are flagged, indicating GPS jumps or data gaps:

$$\text{Flag if: } d_i > 5 \cdot \tilde{d}$$

where $\tilde{d}$ is the median inter-point distance.

#### 1.6 Smooth Velocity Change

Adapted from the Yield Editor SMV (Smooth Velocity) filter. Flags points where the relative speed change between consecutive readings exceeds a user-specified factor:

$$\text{Flag if: } \frac{|v_i - v_{i-1}|}{v_i} > f_{\text{smooth}}$$

*Reference: Sudduth & Drummond (2007), SMV filter.*

#### 1.7 Start/End Delay

**Applicable primarily to yield monitors and other sensors with measurement lag.** Detects ramp-up (filling) and ramp-down (emptying) patterns at pass boundaries. For each transect, the first and last 20 points are scanned. Points below 70% of the transect median value are flagged:

$$\tau_t = 0.70 \cdot \tilde{x}_t$$

$$\text{Flag start points } i \text{ where } x_i < \tau_t \text{ before the first } x_i \geq \tau_t$$

Auto-detection samples 30 transects and enables the filter if >50% show a ramp-up signature (mean of first 5 points < 1.2x mean of next 5 points).

*References:*
> Blackmore, S. (1999). Remedial correction of yield map data. *Precision Agriculture*, *1*(1), 53--66. https://doi.org/10.1023/A:1009969601387

> Sudduth, K. A., & Drummond, S. T. (2007). Yield Editor: Software for removing errors from crop yield maps. *Agronomy Journal*, *99*(6), 1471--1482. https://doi.org/10.2134/agronj2006.0326

---

### Phase 2: Global Outlier Filters

#### 2.1 Adjusted Boxplot using Medcouple (AOI)

Global value bounds are computed using the adjusted boxplot for skewed distributions. Unlike the classical boxplot, which assumes symmetry and can flag >8% of valid observations in skewed data, the adjusted boxplot uses the **medcouple** (MC) to asymmetrically extend the fences.

**Medcouple** (robust measure of skewness, breakdown value = 25%):

$$\mathrm{MC}(X) = \operatorname{med}_{x_i \leq m_X \leq x_j}\; h(x_i, x_j)$$

where:

$$h(x_i, x_j) = \frac{(x_j - m_X) + (x_i - m_X)}{x_j - x_i}$$

**Adjusted fences:**

For MC >= 0 (right-skewed or symmetric):

$$\text{Lower} = Q_1 - k \cdot e^{-4 \cdot \mathrm{MC}} \cdot \mathrm{IQR} \qquad \text{Upper} = Q_3 + k \cdot e^{3 \cdot \mathrm{MC}} \cdot \mathrm{IQR}$$

For MC < 0 (left-skewed):

$$\text{Lower} = Q_1 - k \cdot e^{-3 \cdot \mathrm{MC}} \cdot \mathrm{IQR} \qquad \text{Upper} = Q_3 + k \cdot e^{4 \cdot \mathrm{MC}} \cdot \mathrm{IQR}$$

where $k = 1.5$ (default) and $\mathrm{IQR} = Q_3 - Q_1$. An optional user-adjustable margin widens the bounds further.

*References:*
> Hubert, M., & Vandervieren, E. (2008). An adjusted boxplot for skewed distributions. *Computational Statistics & Data Analysis*, *52*(12), 5186--5201. https://doi.org/10.1016/j.csda.2007.11.008

> Brys, G., Hubert, M., & Struyf, A. (2004). A robust measure of skewness. *Journal of Computational and Graphical Statistics*, *13*(4), 996--1017. https://doi.org/10.1198/106186004X12632

---

### Phase 3: Local / Spatial Filters

#### 3.1 Grid-Based Local Standard Deviation

Points are assigned to grid cells sized as a multiple of the inferred swath width (default: 3x swath). Cell-level mean and standard deviation are computed from **currently clean data only** (points that passed Phases 1 and 2). A point is flagged if it deviates more than *n* standard deviations from its cell mean:

$$\text{Flag if: } |x_p - \bar{x}_C| > n \cdot s_C$$

where $\bar{x}_C$ and $s_C$ are the mean and standard deviation of clean points in cell $C$, and $n$ is the sigma multiplier (default 3.0).

*References:*
> Sudduth, K. A., & Drummond, S. T. (2007). Yield Editor: Software for removing errors from crop yield maps. *Agronomy Journal*, *99*(6), 1471--1482. https://doi.org/10.2134/agronj2006.0326

> Vega, A., Cordoba, M., Castro-Franco, M., & Balzarini, M. (2019). Protocol for automating error removal from yield maps. *Precision Agriculture*, *20*(5), 1030--1044. https://doi.org/10.1007/s11119-018-09632-8

#### 3.2 Cross-Track (Stripy Transect) Detection

**Applicable primarily to yield data and sensors prone to calibration drift.** Compares each transect's median value to a rolling median of 5 neighboring transects. Transects whose ratio deviates beyond thresholds are flagged or corrected:

$$r_t = \frac{\tilde{x}_t}{\tilde{x}_t^{(\mathrm{roll})}}$$

$$\text{Flag if: } r_t > 1.5 \quad \text{or} \quad r_t < 0.66$$

Optional correction scales each transect's values by the inverse ratio:

$$x_i^{(\mathrm{corrected})} = x_i \cdot \frac{\tilde{x}_t^{(\mathrm{roll})}}{\tilde{x}_t}$$

*References:*
> Sudduth, K. A., & Drummond, S. T. (2007). Yield Editor: Software for removing errors from crop yield maps. *Agronomy Journal*, *99*(6), 1471--1482. https://doi.org/10.2134/agronj2006.0326

> Blackmore, S. (1999). Remedial correction of yield map data. *Precision Agriculture*, *1*(1), 53--66. https://doi.org/10.1023/A:1009969601387

---

### Swath Width Inference

Swath width is used to scale grid cells for overlap detection and local filtering. RETA infers it in priority order:

1. From a `swath_width` or `Width` column in the data (if present)
2. From the ratio of convex hull area to total travel distance: $w = A_{\text{hull}} / D_{\text{total}}$
3. Fallback: 10 m

---

## Output Fields

The output layer includes all original fields plus:

| Field | Type | Description |
|---|---|---|
| `dist_m` | Double | Distance to previous point (meters) |
| `bearing_deg` | Double | Travel bearing (degrees, clockwise from north) |
| `speed_mps` | Double | Speed (meters per second) |
| `transect_id` | Integer | Pass/transect identifier (-1 = turn/short) |
| `is_turn` | Boolean | True if point is in a turn zone |
| `is_short` | Boolean | True if point is in a short segment |
| `error_flag` | Integer | Bitmask encoding all triggered filters |
| `error_reason` | String | Human-readable list of triggered filters |
| `filtering_category` | String | Summary category (Clean / Operational Error / Global Outlier / Local Outlier) |
| `Error Type` | String | Same as `filtering_category` |

---

## Further Reading

- Lyle, G., Bryan, B. A., & Ostendorf, B. (2014). Post-processing methods to eliminate erroneous grain yield measurements: Review and directions for future development. *Precision Agriculture*, *15*(4), 377--402. https://doi.org/10.1007/s11119-013-9336-3

- Spekken, M., Anselmi, A. A., & Molin, J. P. (2013). A simple method for filtering spatial data. In J. V. Stafford (Ed.), *Precision Agriculture '13* (pp. 259--266). Wageningen Academic Publishers. https://doi.org/10.3920/978-90-8686-778-3_30

---

## License

This project is licensed under the **GNU General Public License v2.0 or later** (GPL-2.0-or-later). See the [LICENSE](LICENSE) file for the full text.
