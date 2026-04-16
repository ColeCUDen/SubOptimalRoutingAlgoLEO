# routing.py
# GORA and AORA routing algorithms for the LEO satellite network.

# Correlation matrix X indexing convention (follows the paper):
#   rows/cols  0 .. M-1          → gateways   (M = 46)
#   rows/cols  M .. M+N-1        → satellites (N = 720)
#   Total size: (M+N) x (M+N) = 766 x 766

# Edge weight c_{i,j} = current traffic load (Gbps) on the link.
#   - 0.0   : link exists but carries no traffic (valid, low-cost edge)
#   - np.inf: no LOS / link does not exist
# scipy.sparse.csgraph.dijkstra treats 0 as "no edge", so we add EPS to all
# valid edges before calling it and subtract afterwards.


import time
import numpy as np
from scipy.sparse.csgraph import dijkstra as sp_dijkstra
from scipy.sparse import csr_matrix

from .constellation import (
    EARTH_RADIUS_KM,
    build_visibility_matrix,
    gateway_ecef,
)


# This Epsilon offset is necessary so valid zero-load edges survive scipy's
# 0=no-edge convention

_EPS = 1e-9

# Speed of light in km/s
SPEED_OF_LIGHT_KM_S = 299_792.458


# ─── Geometric helpers ────────────────────────────────────────────────────────

# Standard haversine formula that takes two points as (lat, lon) in degrees,
# and returns the great-circle distance in km between them along Earth's surface.
# Used to decide which gateways to include in the sub-graph
def haversine_distance_km(lat1: float, lon1: float,
                          lat2: float, lon2: float) -> float:
    R = EARTH_RADIUS_KM
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dPhi  = phi2 - phi1
    dLam  = np.radians(lon2 - lon1)
    a = np.sin(dPhi / 2)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dLam / 2)**2
    return 2 * R * np.arcsin(np.sqrt(a))


# Finds the point exactly halfway between two surface points. Method: average unit
# Cartesian vectors, normalise, convert back to lat/lon. AORA uses this midpoint
# as the center of its search region. If the two points are antipodal, it returns the
# first point as fallback.
def spherical_midpoint(lat1: float, lon1: float,
                       lat2: float, lon2: float):
    phi1, lam1 = np.radians(lat1), np.radians(lon1)
    phi2, lam2 = np.radians(lat2), np.radians(lon2)

    v1 = np.array([np.cos(phi1)*np.cos(lam1),
                   np.cos(phi1)*np.sin(lam1),
                   np.sin(phi1)])
    v2 = np.array([np.cos(phi2)*np.cos(lam2),
                   np.cos(phi2)*np.sin(lam2),
                   np.sin(phi2)])

    vecMid = v1 + v2
    normMid = np.linalg.norm(vecMid)
    if normMid < 1e-12:
        # antipodal points — return first point as fallback
        return lat1, lon1
    vecMid /= normMid

    latMid = np.degrees(np.arcsin(np.clip(vecMid[2], -1.0, 1.0)))
    lonMid = np.degrees(np.arctan2(vecMid[1], vecMid[0]))
    return latMid, lonMid


# ─── Correlation matrix builder ───────────────────────────────────────────────

# This builds the 766 x 766 matrix that represents the entire network as a graph.
# Indices 0-45 = gateways, indices 46-765 = satellites
# X[i][j] = load in Gbps if gateway i can see satellite j (has line-of-sight)
# X[i][j] = inf if no link exists
# Diagonal is 0 (node to itself)

def build_correlation_matrix(
    gateway_locations,   # list of (lat, lon), length M
    satEcefArr,          # (N, 3) ndarray
    visible,             # (M, N) bool
    rangesKm,            # (M, N) float
    linkLoads=None,      # (M, N) float or None → all zeros
):
    M = len(gateway_locations)
    N = len(satEcefArr)
    size = M + N

    X = np.full((size, size), np.inf)
    np.fill_diagonal(X, 0.0)

    if linkLoads is None:
        linkLoads = np.zeros((M, N), dtype=float)

    # Gateway ↔ Satellite edges — vectorised
    gwIdx, satIdx = np.where(visible)       # indices of all visible (gw, sat) pairs
    vals = linkLoads[gwIdx, satIdx]
    X[gwIdx,       M + satIdx] = vals       # gw row → sat col
    X[M + satIdx,  gwIdx     ] = vals       # sat row → gw col  (symmetric)

    return X


# ─── Dijkstra wrapper ─────────────────────────────────────────────────────────

# Shared Dijkstra wrapper that both GORA and AORA call.
# np.inf entries are treated as absent edges, and valid
# 0.0 entries get +EPS offset so scipy does not ignore them.
# Steps:
# 1. Copy the weight matrix
# 2. Set all inf entries to 0 (scipy's "no edge" convention)
# 3. Add _EPS to all valid edges so zero-load links aren't mistaken for missing edges
# 4. Convert to sparse format and run scipy.sparse.csgraph.dijkstra
# 5. If the target is unreachable, return empty path
# 6. Otherwise, walk backwards through the predecessor array to reconstruct the path
# (target → ... → source, then reverse)
# 7. Subtract the accumulated epsilon from the total cost to get the true cost
# Returns: (path as list of node indices, total cost)

def dijkstra_on_matrix(weightMatrix: np.ndarray, sourceIdx: int, targetIdx: int):
    n = weightMatrix.shape[0]

    # Build scipy-compatible matrix: 0 means no edge
    W = weightMatrix.copy()
    finiteMask = np.isfinite(W) & (W >= 0)
    W[~finiteMask] = 0.0                  # absent edges → 0 (scipy convention)
    W[finiteMask] += _EPS                 # valid edges  → >0

    np.fill_diagonal(W, 0.0)             # self-loops stay 0 (ignored by scipy)

    sparseW = csr_matrix(W)
    distArr, predArr = sp_dijkstra(
        sparseW, directed=False,
        indices=sourceIdx,
        return_predecessors=True,
    )

    if np.isinf(distArr[targetIdx]):
        return [], np.inf

    # Reconstruct path
    path = []
    cur = targetIdx
    while cur != sourceIdx:
        path.append(cur)
        cur = predArr[cur]
        if cur < 0:
            return [], np.inf
    path.append(sourceIdx)
    path.reverse()

    # Total cost = sum of original (non-EPS) edge weights
    totalCost = distArr[targetIdx] - _EPS * (len(path) - 1)
    return path, max(totalCost, 0.0)


# Takes a path (list of node indices) and computes the real-world signal delay in milliseconds.
# For each consecutive pair of nodes in the path, it:
# 1. Looks up their ECEF (x,y,z) coordinates — gateway positions if index < 46, satellite positions otherwise
# 2. Computes the Euclidean distance between them (slant range)
# 3. Divides by speed of light, and multiplies by 1000 to get ms
def path_propagation_delay_ms(path: list, weightMatrix: np.ndarray,
                               gatewayEcefArr: np.ndarray,
                               satEcefArr: np.ndarray,
                               numGateways: int) -> float:
    if len(path) < 2:
        return 0.0

    def nodeEcef(idx):
        if idx < numGateways:
            return gatewayEcefArr[idx]
        return satEcefArr[idx - numGateways]

    totalMs = 0.0
    for a, b in zip(path[:-1], path[1:]):
        distKm = np.linalg.norm(nodeEcef(b) - nodeEcef(a))
        totalMs += distKm / SPEED_OF_LIGHT_KM_S * 1000.0
    return totalMs


# ─── GORA ─────────────────────────────────────────────────────────────────────

# Global Optimal Routing Algorithm. Just computes Dijkstra's on the whole matrix.
# Complexity: O((M+N)^2) per query.
# route(srcGw, dstGw) starts a timer, calls dijkstra_on_matrix on the full matrix,
# stops the timer, computes the propagation delay on the resulting path, and returns
# the path, cost, matrixSize, delayMs, and compTimeMs.
class GORA:

    def __init__(self, correlationMatrix: np.ndarray,
                 gatewayEcefArr: np.ndarray,
                 satEcefArr: np.ndarray,
                 numGateways: int):
        self.X           = correlationMatrix
        self.gwEcef      = gatewayEcefArr
        self.satEcef     = satEcefArr
        self.M           = numGateways
        self.matrixSize  = correlationMatrix.shape[0]   # M + N

    def route(self, srcGw: int, dstGw: int):
        t0 = time.perf_counter()
        path, cost = dijkstra_on_matrix(self.X, srcGw, dstGw)
        compTimeMs = (time.perf_counter() - t0) * 1000.0

        delayMs = path_propagation_delay_ms(
            path, self.X, self.gwEcef, self.satEcef, self.M
        ) if path else np.inf

        return path, cost, self.matrixSize, delayMs, compTimeMs


# ─── AORA ─────────────────────────────────────────────────────────────────────

# Approximate Optimal Routing Algorithm.
# Reduces the correlation matrix by selecting only gateways near the
# spherical midpoint of the source–destination pair, then runs Dijkstra
# on the smaller sub-matrix.  Falls back to GORA if no path is found.
#
# _selectSubgraphIndices(srcGw, dstGw) selects the gateway and satellite
# indices for the sub-correlation matrix.
# Steps:
# 1. Compute spherical midpoint k of src and dst.
# 2. dKp = haversine(k, srcGw).
# 3. Radius R = dKp + extraRadiusKm.
# 4. Include all gateways g where haversine(k, g) <= R. (src and dst are always included.)
# 5. Include all satellites visible from any selected gateway.
# Returns: gwIndices  : sorted list of gateway indices in the sub-matrix
#          satIndices : sorted list of satellite indices in the sub-matrix
#
# _buildSubmatrix(gwIndices, satIndices) extracts the sub-correlation matrix and returns it
# and the index mappings for the gateways and satellites.
#
# route(srcGw, dstGw) routes between two gateway indices using the sub-matrix Dijkstra. It
# falls back to full matrix (GORA) if there is no path found in the sub-matrix. It returns:
# path, cost, subSize, delayMs, compTimeMs, usedFallback

class AORA:

    def __init__(self, correlationMatrix: np.ndarray,
                 gatewayLocations,              # list of (lat, lon)
                 gatewayEcefArr: np.ndarray,
                 satEcefArr: np.ndarray,
                 visible: np.ndarray,
                 numGateways: int,
                 extraRadiusKm: float = 600.0):
        self.X            = correlationMatrix
        self.gwLocs       = gatewayLocations
        self.gwEcef       = gatewayEcefArr
        self.satEcef      = satEcefArr
        self.visible      = visible           # (M, N) bool
        self.M            = numGateways
        self.N            = len(satEcefArr)
        self.extraRadius  = extraRadiusKm

    def _selectSubgraphIndices(self, srcGw: int, dstGw: int):
        latS, lonS = self.gwLocs[srcGw]
        latD, lonD = self.gwLocs[dstGw]

        latK, lonK = spherical_midpoint(latS, lonS, latD, lonD)
        dKp = haversine_distance_km(latK, lonK, latS, lonS)
        R   = dKp + self.extraRadius

        # Select gateways within radius R of midpoint
        gwSet = set([srcGw, dstGw])
        for g, (latG, lonG) in enumerate(self.gwLocs):
            if haversine_distance_km(latK, lonK, latG, lonG) <= R:
                gwSet.add(g)
        gwIndices = sorted(gwSet)

        # Select satellites visible from any selected gateway
        satSet = set()
        for g in gwIndices:
            satSet.update(np.where(self.visible[g])[0].tolist())
        satIndices = sorted(satSet)

        return gwIndices, satIndices

    def _buildSubmatrix(self, gwIndices, satIndices):
        m = len(gwIndices)
        n = len(satIndices)
        size = m + n
        subX = np.full((size, size), np.inf)
        np.fill_diagonal(subX, 0.0)

        gwMap  = {g: i       for i, g in enumerate(gwIndices)}
        satMap = {s: m + j   for j, s in enumerate(satIndices)}

        for g in gwIndices:
            for s in satIndices:
                if self.visible[g, s]:
                    localGw  = gwMap[g]
                    localSat = satMap[s]
                    val = self.X[g, self.M + s]
                    subX[localGw, localSat] = val
                    subX[localSat, localGw] = val

        return subX, gwMap, satMap

    def route(self, srcGw: int, dstGw: int):
        t0 = time.perf_counter()

        gwIndices, satIndices = self._selectSubgraphIndices(srcGw, dstGw)
        subX, gwMap, satMap   = self._buildSubmatrix(gwIndices, satIndices)

        subSrc  = gwMap[srcGw]
        subDst  = gwMap[dstGw]
        subSize = subX.shape[0]

        subPath, cost = dijkstra_on_matrix(subX, subSrc, subDst)

        usedFallback = False

        if not subPath:
            # Fallback: use full correlation matrix
            usedFallback = True
            subSize = self.X.shape[0]
            globalPath, cost = dijkstra_on_matrix(self.X, srcGw, dstGw)
        else:
            # Map local indices back to global indices
            invGw  = {v: k for k, v in gwMap.items()}
            invSat = {v: k for k, v in satMap.items()}

            globalPath = []
            for localIdx in subPath:
                if localIdx in invGw:
                    globalPath.append(invGw[localIdx])
                else:
                    # satellite: local → global satellite idx → global matrix idx
                    globalSat = invSat[localIdx]
                    globalPath.append(self.M + globalSat)

        compTimeMs = (time.perf_counter() - t0) * 1000.0

        delayMs = path_propagation_delay_ms(
            globalPath, self.X, self.gwEcef, self.satEcef, self.M
        ) if globalPath else np.inf

        return globalPath, cost, subSize, delayMs, compTimeMs, usedFallback
