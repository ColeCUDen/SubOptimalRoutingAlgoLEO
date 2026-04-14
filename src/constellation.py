import numpy as np

# Orbital Parameters
EARTH_RADIUS_KM = 6371.0
MU_KM3_S2 = 3.986004418e5 # Earth gravitational parameter

NUM_PLANES = 18
SATS_PER_PLANE = 40
NUM_SATELLITES  = NUM_PLANES * SATS_PER_PLANE # 720

ALTITUDE_KM = 1200.0
INCLINATION_DEG = 90.0 # polar orbits
PHASING_DEG = 4.5 # mean-anomaly offset between adjacent planes
MIN_ELEVATION_DEG = 10.0 # link visibility threshold

ORBIT_RADIUS_KM = EARTH_RADIUS_KM + ALTITUDE_KM  # 7571 km
ORBITAL_PERIOD_S = 2 * np.pi * np.sqrt(ORBIT_RADIUS_KM**3 / MU_KM3_S2) # ~6244 s
MEAN_MOTION_RAD_S = 2 * np.pi / ORBITAL_PERIOD_S


def satellite_position_ecef(plane_idx: int, sat_idx: int, time_s: float = 0.0):
    inc  = np.radians(INCLINATION_DEG)
    raan = np.radians(plane_idx * 360.0 / NUM_PLANES)

    M0 = np.radians(sat_idx  * 360.0 / SATS_PER_PLANE + plane_idx * PHASING_DEG)
    u = M0 + MEAN_MOTION_RAD_S * time_s   # current argument of latitude

    # Position in orbital plane
    x_orb = ORBIT_RADIUS_KM * np.cos(u)
    y_orb = ORBIT_RADIUS_KM * np.sin(u)

    # Rotate: first by inclination (tilt orbital plane), then by RAAN
    # Rotation matrix R_z(RAAN) * R_x(inc)
    ci, si = np.cos(inc),  np.sin(inc)
    cr, sr = np.cos(raan), np.sin(raan)

    x = cr * x_orb - sr * ci * y_orb
    y = sr * x_orb + cr * ci * y_orb
    z = si * y_orb

    return np.array([x, y, z])


def gateway_ecef(lat_deg: float, lon_deg: float) -> np.ndarray:
    # Convert geodetic (lat, lon) on Earth's surface to ECEF km
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    x = EARTH_RADIUS_KM * np.cos(lat) * np.cos(lon)
    y = EARTH_RADIUS_KM * np.cos(lat) * np.sin(lon)
    z = EARTH_RADIUS_KM * np.sin(lat)
    return np.array([x, y, z])


def elevation_angle_deg(gw_ecef: np.ndarray, sat_ecef: np.ndarray) -> float:
    los = sat_ecef - gw_ecef
    los_norm = los / np.linalg.norm(los)
    up  = gw_ecef / np.linalg.norm(gw_ecef)
    sin_el = np.dot(los_norm, up)
    return np.degrees(np.arcsin(np.clip(sin_el, -1.0, 1.0)))


def slant_range_km(gw_ecef: np.ndarray, sat_ecef: np.ndarray) -> float:
    # Slant range (km) between gateway and satellite
    return float(np.linalg.norm(sat_ecef - gw_ecef))


def build_all_satellite_positions(time_s: float = 0.0) -> np.ndarray:
    positions = np.empty((NUM_SATELLITES, 3))
    for p in range(NUM_PLANES):
        for s in range(SATS_PER_PLANE):
            idx = p * SATS_PER_PLANE + s
            positions[idx] = satellite_position_ecef(p, s, time_s)
    return positions


def build_visibility_matrix(gateway_ecef_list: np.ndarray, satellite_ecef_arr: np.ndarray,min_el: float = MIN_ELEVATION_DEG):
    M = len(gateway_ecef_list)
    N = len(satellite_ecef_arr)
    visible = np.zeros((M, N), dtype=bool)
    ranges  = np.zeros((M, N), dtype=float)

    for i, gw in enumerate(gateway_ecef_list):
        los_vecs = satellite_ecef_arr - gw # (N, 3)
        dists = np.linalg.norm(los_vecs, axis=1) # (N,)
        up = gw / np.linalg.norm(gw)
        sin_el = (los_vecs / dists[:, None]) @ up # (N,) dot products
        el_deg = np.degrees(np.arcsin(np.clip(sin_el, -1.0, 1.0)))
        mask = el_deg >= min_el
        visible[i] = mask
        ranges[i][mask] = dists[mask]

    return visible, ranges