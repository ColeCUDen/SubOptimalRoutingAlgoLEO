"""
routing.py
GORA and AORA routing algorithms for the LEO satellite network.

Correlation matrix X indexing convention (follows the paper):
  rows/cols  0 .. M-1          → gateways   (M = 46)
  rows/cols  M .. M+N-1        → satellites (N = 720)
  Total size: (M+N) x (M+N) = 766 x 766

Edge weight c_{i,j} = current traffic load (Gbps) on the link.
  - 0.0   : link exists but carries no traffic (valid, low-cost edge)
  - np.inf: no LOS / link does not exist
scipy.sparse.csgraph.dijkstra treats 0 as "no edge", so we add EPS to all
valid edges before calling it and subtract afterwards.
"""

import time
import numpy as np
from scipy.sparse.csgraph import dijkstra as sp_dijkstra
from scipy.sparse import csr_matrix

from .constellation import (
    EARTH_RADIUS_KM,
    build_visibility_matrix,
    gateway_ecef,
)

# Epsilon offset so valid zero-load edges survive scipy's 0=no-edge convention
_EPS = 1e-9

# Speed of light (km/s)
SPEED_OF_LIGHT_KM_S = 299_792.458


# ─── Geometric helpers ────────────────────────────────────────────────────────

def haversine_distance_km(lat1: float, lon1: float,
                          lat2: float, lon2: float) -> float:
    """Great-circle distance (km) between two (lat, lon) points in degrees."""
    R = EARTH_RADIUS_KM
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi  = phi2 - phi1
    dlam  = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


def spherical_midpoint(lat1: float, lon1: float,
                       lat2: float, lon2: float):
    """
    Spherical midpoint between two surface points.
    Method: average unit Cartesian vectors, normalise, convert back to lat/lon.
    Returns (lat_deg, lon_deg).
    """
    phi1, lam1 = np.radians(lat1), np.radians(lon1)
    phi2, lam2 = np.radians(lat2), np.radians(lon2)

    v1 = np.array([np.cos(phi1)*np.cos(lam1),
                   np.cos(phi1)*np.sin(lam1),
                   np.sin(phi1)])
    v2 = np.array([np.cos(phi2)*np.cos(lam2),
                   np.cos(phi2)*np.sin(lam2),
                   np.sin(phi2)])

    vm = v1 + v2
    nm = np.linalg.norm(vm)
    if nm < 1e-12:
        # antipodal points — return first point as fallback
        return lat1, lon1
    vm /= nm

    lat_m = np.degrees(np.arcsin(np.clip(vm[2], -1.0, 1.0)))
    lon_m = np.degrees(np.arctan2(vm[1], vm[0]))
    return lat_m, lon_m


# ─── Correlation matrix builder ───────────────────────────────────────────────

def build_correlation_matrix(
    gateway_locations,   # list of (lat, lon), length M
    sat_ecef_arr,        # (N, 3) ndarray
    visible,             # (M, N) bool
    ranges_km,           # (M, N) float
    link_loads=None,     # (M, N) float or None → all zeros
):
    """
    Build the (M+N)×(M+N) correlation matrix X.

    X[i, j] = link load (Gbps)  if link exists (LOS)
    X[i, j] = np.inf             if no link
    Diagonal is 0.

    Index mapping:
      gateway  i  →  row/col i          (0 .. M-1)
      satellite s  →  row/col M + s     (M .. M+N-1)
    """
    M = len(gateway_locations)
    N = len(sat_ecef_arr)
    size = M + N

    X = np.full((size, size), np.inf)
    np.fill_diagonal(X, 0.0)

    if link_loads is None:
        link_loads = np.zeros((M, N), dtype=float)

    # Gateway ↔ Satellite edges — vectorised
    gi, si = np.where(visible)           # indices of all visible (gw, sat) pairs
    vals   = link_loads[gi, si]
    X[gi,      M + si] = vals            # gw row → sat col
    X[M + si,  gi    ] = vals            # sat row → gw col  (symmetric)

    return X


# ─── Dijkstra wrapper ─────────────────────────────────────────────────────────

def dijkstra_on_matrix(weight_matrix: np.ndarray, source_idx: int, target_idx: int):
    """
    Run Dijkstra on a dense weight matrix.
    - np.inf entries are treated as absent edges.
    - Valid 0.0 entries get +EPS offset so scipy does not ignore them.

    Returns
    -------
    path       : list of node indices from source to target, or [] if unreachable
    total_cost : float — sum of original edge weights along path
    """
    n = weight_matrix.shape[0]

    # Build scipy-compatible matrix: 0 means no edge
    W = weight_matrix.copy()
    finite_mask = np.isfinite(W) & (W >= 0)
    W[~finite_mask] = 0.0                  # absent edges → 0 (scipy convention)
    W[finite_mask] += _EPS                 # valid edges  → >0

    np.fill_diagonal(W, 0.0)              # self-loops stay 0 (ignored by scipy)

    sparse_W = csr_matrix(W)
    dist_arr, pred_arr = sp_dijkstra(
        sparse_W, directed=False,
        indices=source_idx,
        return_predecessors=True,
    )

    if np.isinf(dist_arr[target_idx]):
        return [], np.inf

    # Reconstruct path
    path = []
    cur = target_idx
    while cur != source_idx:
        path.append(cur)
        cur = pred_arr[cur]
        if cur < 0:
            return [], np.inf
    path.append(source_idx)
    path.reverse()

    # Total cost = sum of original (non-EPS) edge weights
    total_cost = dist_arr[target_idx] - _EPS * (len(path) - 1)
    return path, max(total_cost, 0.0)


def path_propagation_delay_ms(path: list, weight_matrix: np.ndarray,
                               gateway_ecef_arr: np.ndarray,
                               sat_ecef_arr: np.ndarray,
                               num_gateways: int) -> float:
    """
    Compute total propagation delay (ms) along a path.
    Delay per hop = slant_range_km / speed_of_light_km_s * 1000 ms.
    """
    if len(path) < 2:
        return 0.0

    def node_ecef(idx):
        if idx < num_gateways:
            return gateway_ecef_arr[idx]
        return sat_ecef_arr[idx - num_gateways]

    total_ms = 0.0
    for a, b in zip(path[:-1], path[1:]):
        dist_km   = np.linalg.norm(node_ecef(b) - node_ecef(a))
        total_ms += dist_km / SPEED_OF_LIGHT_KM_S * 1000.0
    return total_ms


# ─── GORA ─────────────────────────────────────────────────────────────────────

class GORA:
    """
    Global Optimal Routing Algorithm.
    Runs Dijkstra on the full (M+N)×(M+N) correlation matrix.
    Complexity: O((M+N)^2) per query.
    """

    def __init__(self, correlation_matrix: np.ndarray,
                 gateway_ecef_arr: np.ndarray,
                 sat_ecef_arr: np.ndarray,
                 num_gateways: int):
        self.X             = correlation_matrix
        self.gw_ecef       = gateway_ecef_arr
        self.sat_ecef      = sat_ecef_arr
        self.M             = num_gateways
        self.matrix_size   = correlation_matrix.shape[0]   # M + N

    def route(self, src_gw: int, dst_gw: int):
        """
        Route between two gateway indices.

        Returns
        -------
        path         : list of node indices
        total_cost   : float
        matrix_size  : int  (always M+N = 766)
        delay_ms     : float
        comp_time_ms : float
        """
        t0 = time.perf_counter()
        path, cost = dijkstra_on_matrix(self.X, src_gw, dst_gw)
        comp_time_ms = (time.perf_counter() - t0) * 1000.0

        delay_ms = path_propagation_delay_ms(
            path, self.X, self.gw_ecef, self.sat_ecef, self.M
        ) if path else np.inf

        return path, cost, self.matrix_size, delay_ms, comp_time_ms


# ─── AORA ─────────────────────────────────────────────────────────────────────

class AORA:
    """
    Approximate Optimal Routing Algorithm.
    Reduces the correlation matrix by selecting only gateways near the
    spherical midpoint of the source–destination pair, then runs Dijkstra
    on the smaller sub-matrix.  Falls back to GORA if no path is found.
    """

    def __init__(self, correlation_matrix: np.ndarray,
                 gateway_locations,          # list of (lat, lon)
                 gateway_ecef_arr: np.ndarray,
                 sat_ecef_arr: np.ndarray,
                 visible: np.ndarray,
                 num_gateways: int,
                 extra_radius_km: float = 600.0):
        self.X            = correlation_matrix
        self.gw_locs      = gateway_locations
        self.gw_ecef      = gateway_ecef_arr
        self.sat_ecef     = sat_ecef_arr
        self.visible      = visible           # (M, N) bool
        self.M            = num_gateways
        self.N            = len(sat_ecef_arr)
        self.extra_radius = extra_radius_km

    def _select_subgraph_indices(self, src_gw: int, dst_gw: int):
        """
        Select gateway and satellite indices for the sub-correlation matrix.

        Steps (from the paper):
          1. Compute spherical midpoint k of src and dst.
          2. D_kp = haversine(k, src_gw).
          3. Radius R = D_kp + extra_radius_km.
          4. Include all gateways g where haversine(k, g) <= R.
             (src and dst are always included.)
          5. Include all satellites visible from any selected gateway.

        Returns
        -------
        gw_indices  : sorted list of gateway indices in the sub-matrix
        sat_indices : sorted list of satellite indices in the sub-matrix
        """
        lat_s, lon_s = self.gw_locs[src_gw]
        lat_d, lon_d = self.gw_locs[dst_gw]

        lat_k, lon_k = spherical_midpoint(lat_s, lon_s, lat_d, lon_d)
        D_kp = haversine_distance_km(lat_k, lon_k, lat_s, lon_s)
        R    = D_kp + self.extra_radius

        # Select gateways within radius R of midpoint
        gw_set = set([src_gw, dst_gw])
        for g, (lat_g, lon_g) in enumerate(self.gw_locs):
            if haversine_distance_km(lat_k, lon_k, lat_g, lon_g) <= R:
                gw_set.add(g)
        gw_indices = sorted(gw_set)

        # Select satellites visible from any selected gateway
        sat_set = set()
        for g in gw_indices:
            sat_set.update(np.where(self.visible[g])[0].tolist())
        sat_indices = sorted(sat_set)

        return gw_indices, sat_indices

    def _build_submatrix(self, gw_indices, sat_indices):
        """
        Extract the sub-correlation matrix and return index mappings.

        Returns
        -------
        sub_X      : (m+n) × (m+n) sub-matrix
        gw_map     : dict  global_gw_idx  → local_idx
        sat_map    : dict  global_sat_idx → local_idx (offset by m)
        """
        m = len(gw_indices)
        n = len(sat_indices)
        size = m + n
        sub_X = np.full((size, size), np.inf)
        np.fill_diagonal(sub_X, 0.0)

        gw_map  = {g: i       for i, g in enumerate(gw_indices)}
        sat_map = {s: m + j   for j, s in enumerate(sat_indices)}

        for g in gw_indices:
            for s in sat_indices:
                if self.visible[g, s]:
                    li = gw_map[g]
                    lj = sat_map[s]
                    val = self.X[g, self.M + s]
                    sub_X[li, lj] = val
                    sub_X[lj, li] = val

        return sub_X, gw_map, sat_map

    def route(self, src_gw: int, dst_gw: int):
        """
        Route between two gateway indices using sub-matrix Dijkstra.
        Falls back to full matrix (GORA) if no path found in sub-matrix.

        Returns
        -------
        path          : list of *global* node indices
        total_cost    : float
        sub_size      : int   (size of sub-matrix used)
        delay_ms      : float
        comp_time_ms  : float
        used_fallback : bool
        """
        t0 = time.perf_counter()

        gw_indices, sat_indices = self._select_subgraph_indices(src_gw, dst_gw)
        sub_X, gw_map, sat_map  = self._build_submatrix(gw_indices, sat_indices)

        sub_src = gw_map[src_gw]
        sub_dst = gw_map[dst_gw]
        sub_size = sub_X.shape[0]

        sub_path, cost = dijkstra_on_matrix(sub_X, sub_src, sub_dst)

        used_fallback = False

        if not sub_path:
            # Fallback: use full correlation matrix
            used_fallback = True
            sub_size = self.X.shape[0]
            global_path, cost = dijkstra_on_matrix(self.X, src_gw, dst_gw)
        else:
            # Map local indices back to global indices
            inv_gw  = {v: k for k, v in gw_map.items()}
            inv_sat = {v: k for k, v in sat_map.items()}

            global_path = []
            for local_idx in sub_path:
                if local_idx in inv_gw:
                    global_path.append(inv_gw[local_idx])
                else:
                    # satellite: local → global satellite idx → global matrix idx
                    global_sat = inv_sat[local_idx]
                    global_path.append(self.M + global_sat)

        comp_time_ms = (time.perf_counter() - t0) * 1000.0

        delay_ms = path_propagation_delay_ms(
            global_path, self.X, self.gw_ecef, self.sat_ecef, self.M
        ) if global_path else np.inf

        return global_path, cost, sub_size, delay_ms, comp_time_ms, used_fallback
