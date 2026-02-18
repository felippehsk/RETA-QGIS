"""
CRS-aware distance, area, and coordinate conversion. No preset scale numbers;
uses pyproj for proper projection/ellipsoid calculations. Works with CRS authid
strings (e.g. EPSG:4326, EPSG:32615). When running in QGIS, pass source.sourceCrs().authid().
"""

import numpy as np

try:
    from pyproj import CRS, Geod, Transformer
    PYPROJ_AVAILABLE = True
except ImportError:
    PYPROJ_AVAILABLE = False


def _normalize_crs(crs_spec):
    """QGIS QgsCoordinateReferenceSystem or authid string -> authid string."""
    if crs_spec is None:
        return None
    if hasattr(crs_spec, "authid"):
        try:
            return crs_spec.authid() if getattr(crs_spec, "isValid", lambda: True)() else None
        except Exception:
            return None
    return str(crs_spec) if crs_spec else None


def _crs_from_user_input(crs_spec):
    """crs_spec: authid string (e.g. EPSG:4326), or QGIS CRS, or None."""
    spec = _normalize_crs(crs_spec)
    if spec is None or not PYPROJ_AVAILABLE:
        return None
    try:
        return CRS.from_user_input(spec)
    except Exception:
        return None


def is_geographic(crs_spec):
    """True if CRS is geographic (deg), False if projected (e.g. UTM in m)."""
    spec = _normalize_crs(crs_spec)
    crs = _crs_from_user_input(spec)
    if crs is None:
        return True  # default assume WGS84 degrees
    try:
        return crs.is_geographic
    except Exception:
        return True


def _utm_zone_epsg(lon_deg, lat_deg):
    """Approximate UTM zone EPSG (northern hemisphere). lon, lat in degrees."""
    zone = int((float(lon_deg) + 180) / 6) + 1
    zone = max(1, min(60, zone))
    # Northern hemisphere
    return 32600 + zone if float(lat_deg) >= 0 else 32700 + zone


def distances_m(lon1, lat1, lon2, lat2, crs_spec):
    """
    Vectorized distance in meters between (lon1, lat1) and (lon2, lat2).
    lon1, lat1, lon2, lat2: arrays or scalars. Returns same shape.
    """
    if crs_spec is None or not PYPROJ_AVAILABLE:
        # Fallback: Haversine (WGS84)
        lon1, lat1, lon2, lat2 = np.atleast_1d(lon1, lat1, lon2, lat2)
        R = 6371000.0
        phi1 = np.radians(lat1)
        phi2 = np.radians(lat2)
        dphi = np.radians(lat2 - lat1)
        dlam = np.radians(lon2 - lon1)
        a = np.sin(dphi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam / 2) ** 2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        return (R * c).reshape(np.broadcast_shapes(lon1.shape, lat1.shape, lon2.shape, lat2.shape))

    spec = _normalize_crs(crs_spec)
    crs = _crs_from_user_input(spec)
    if crs is None:
        return distances_m(lon1, lat1, lon2, lat2, None)

    if crs.is_geographic:
        geod = Geod(ellps="WGS84")
        lon1, lat1, lon2, lat2 = np.atleast_1d(lon1, lat1, lon2, lat2)
        _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
        return dist.reshape(np.broadcast_shapes(lon1.shape, lat1.shape, lon2.shape, lat2.shape))
    else:
        # Projected: units are usually meters
        dx = np.atleast_1d(lon2 - lon1).astype(float)
        dy = np.atleast_1d(lat2 - lat1).astype(float)
        return np.sqrt(dx * dx + dy * dy).reshape(np.broadcast_shapes(dx.shape, dy.shape))


def points_to_meters(lons, lats, crs_spec):
    """
    Convert point arrays to meters (x_m, y_m) for grid/area in m.
    If geographic: transform to local UTM. If projected: return as-is (assume meters).
    """
    spec = _normalize_crs(crs_spec)
    lons = np.atleast_1d(np.asarray(lons, dtype=float))
    lats = np.atleast_1d(np.asarray(lats, dtype=float))
    if spec is None or not PYPROJ_AVAILABLE:
        # Approximate: degrees to meters at centroid (no preset scale in API; internal fallback)
        lon0, lat0 = lons.mean(), lats.mean()
        R = 6371000.0
        scale_y = np.pi * R / 180.0
        scale_x = scale_y * np.cos(np.radians(lat0))
        x_m = (lons - lon0) * scale_x
        y_m = (lats - lat0) * scale_y
        return x_m, y_m

    crs = _crs_from_user_input(spec)
    if crs is None:
        return points_to_meters(lons, lats, None)

    if crs.is_geographic:
        zone_epsg = _utm_zone_epsg(lons.mean(), lats.mean())
        try:
            trans = Transformer.from_crs(spec, f"EPSG:{zone_epsg}", always_xy=True)
            x_m, y_m = trans.transform(lons, lats)
            return np.asarray(x_m, dtype=float), np.asarray(y_m, dtype=float)
        except Exception:
            # Fallback
            return points_to_meters(lons, lats, None)
    else:
        return np.asarray(lons, dtype=float), np.asarray(lats, dtype=float)


def area_m2_from_points(lons, lats, crs_spec):
    """
    Area in m^2 from point cloud: convex hull if available, else bounding box.
    Uses points_to_meters so no preset scale; CRS-aware.
    """
    x_m, y_m = points_to_meters(lons, lats, crs_spec)
    if len(x_m) < 2:
        return 0.0
    try:
        from scipy.spatial import ConvexHull
        pts = np.column_stack([x_m, y_m])
        hull = ConvexHull(pts)
        idx = np.append(hull.vertices, hull.vertices[0])
        return 0.5 * abs(np.sum(x_m[idx[:-1]] * y_m[idx[1:]] - x_m[idx[1:]] * y_m[idx[:-1]]))
    except Exception:
        return float(np.ptp(x_m) * np.ptp(y_m))
