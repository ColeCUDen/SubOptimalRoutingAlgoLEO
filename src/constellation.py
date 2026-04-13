"""
constellation.py
LEO satellite constellation parameters and position/visibility computations.
Constellation: 720 satellites, 18 polar planes x 40 sats/plane, 1200 km altitude.
"""

import numpy as np

# ─── Orbital Parameters ───────────────────────────────────────────────────────
EARTH_RADIUS_KM   = 6371.0
MU_KM3_S2         = 3.986004418e5   # Earth gravitational parameter

NUM_PLANES        = 18
SATS_PER_PLANE    = 40
NUM_SATELLITES    = NUM_PLANES * SATS_PER_PLANE   # 720

ALTITUDE_KM       = 1200.0
INCLINATION_DEG   = 90.0            # polar orbits
PHASING_DEG       = 4.5             # mean-anomaly offset between adjacent planes
MIN_ELEVATION_DEG = 10.0            # link visibility threshold

ORBIT_RADIUS_KM   = EARTH_RADIUS_KM + ALTITUDE_KM  # 7571 km
ORBITAL_PERIOD_S  = 2 * np.pi * np.sqrt(ORBIT_RADIUS_KM**3 / MU_KM3_S2)  # ~6244 s
MEAN_MOTION_RAD_S = 2 * np.pi / ORBITAL_PERIOD_S


def satellite_position_ecef(plane_idx: int, sat_idx: int, time_s: float = 0.0):
    """
    Return ECEF Cartesian coordinates (x, y, z) in km for one satellite.

    Orbital model (analytical circular orbit):
      - RAAN  = plane_idx * (360 / NUM_PLANES)  degrees
      - M0    = sat_idx  * (360 / SATS_PER_PLANE)
                + plane_idx * PHASING_DEG        degrees  (Walker phasing)
      - u(t)  = M0 + n*t  (mean anomaly = true anomaly for circular orbit)
      - Satellite position in orbital plane (perifocal):
            r_orb = [R*cos(u), R*sin(u), 0]
      - Rotate by inclination (i) around x-axis, then RAAN (Ω) around z-axis.
    """
    inc  = np.radians(INCLINATION_DEG)
    raan = np.radians(plane_idx * 360.0 / NUM_PLANES)

    M0   = np.radians(sat_idx  * 360.0 / SATS_PER_PLANE
                      + plane_idx * PHASING_DEG)
    u    = M0 + MEAN_MOTION_RAD_S * time_s   # current argument of latitude

    # Position in orbital plane
    x_orb = ORBIT_RADIUS_KM * np.cos(u)
    y_orb = ORBIT_RADIUS_KM * np.sin(u)

    # Rotate: first by inclination (tilt orbital plane), then by RAAN
    # Rotation matrix R_z(RAAN) * R_x(inc)
    ci, si = np.cos(inc),  np.sin(inc)
    cr, sr = np.cos(raan), np.sin(raan)

    x = cr * x_orb - sr * ci * y_orb
    y = sr * x_orb + cr * ci * y_orb
    z =              si * y_orb

    return np.array([x, y, z])


def gateway_ecef(lat_deg: float, lon_deg: float) -> np.ndarray:
    """Convert geodetic (lat, lon) on Earth's surface to ECEF km."""
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    x = EARTH_RADIUS_KM * np.cos(lat) * np.cos(lon)
    y = EARTH_RADIUS_KM * np.cos(lat) * np.sin(lon)
    z = EARTH_RADIUS_KM * np.sin(lat)
    return np.array([x, y, z])


def elevation_angle_deg(gw_ecef: np.ndarray, sat_ecef: np.ndarray) -> float:
    """
    Elevation angle (degrees) of satellite as seen from a ground gateway.
    el = arcsin( dot(LOS_unit, up_unit) )
    where up_unit = gw_ecef / |gw_ecef|  (outward normal at gateway).
    """
    los = sat_ecef - gw_ecef
    los_norm = los / np.linalg.norm(los)
    up  = gw_ecef / np.linalg.norm(gw_ecef)
    sin_el = np.dot(los_norm, up)
    return np.degrees(np.arcsin(np.clip(sin_el, -1.0, 1.0)))


def slant_range_km(gw_ecef: np.ndarray, sat_ecef: np.ndarray) -> float:
    """Slant range (km) between gateway and satellite."""
    return float(np.linalg.norm(sat_ecef - gw_ecef))


def build_all_satellite_positions(time_s: float = 0.0) -> np.ndarray:
    """
    Compute ECEF positions for all 720 satellites.

    Returns
    -------
    positions : ndarray of shape (720, 3)
        Row i = satellite index i = plane_idx*SATS_PER_PLANE + sat_idx
    """
    positions = np.empty((NUM_SATELLITES, 3))
    for p in range(NUM_PLANES):
        for s in range(SATS_PER_PLANE):
            idx = p * SATS_PER_PLANE + s
            positions[idx] = satellite_position_ecef(p, s, time_s)
    return positions


def build_visibility_matrix(
    gateway_ecef_list: np.ndarray,
    satellite_ecef_arr: np.ndarray,
    min_el: float = MIN_ELEVATION_DEG,
):
    """
    Compute visibility (boolean) and slant-range matrices.

    Parameters
    ----------
    gateway_ecef_list : (M, 3) array  — gateway ECEF positions
    satellite_ecef_arr: (N, 3) array  — satellite ECEF positions
    min_el            : minimum elevation angle in degrees

    Returns
    -------
    visible : (M, N) bool array
    ranges  : (M, N) float array — slant ranges in km (0 where not visible)
    """
    M = len(gateway_ecef_list)
    N = len(satellite_ecef_arr)
    visible = np.zeros((M, N), dtype=bool)
    ranges  = np.zeros((M, N), dtype=float)

    for i, gw in enumerate(gateway_ecef_list):
        los_vecs = satellite_ecef_arr - gw           # (N, 3)
        dists    = np.linalg.norm(los_vecs, axis=1)  # (N,)
        up       = gw / np.linalg.norm(gw)
        sin_el   = (los_vecs / dists[:, None]) @ up  # (N,) dot products
        el_deg   = np.degrees(np.arcsin(np.clip(sin_el, -1.0, 1.0)))
        mask     = el_deg >= min_el
        visible[i]         = mask
        ranges[i][mask]    = dists[mask]

    return visible, ranges
