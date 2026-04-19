import time
import random
import numpy as np

from .constellation import (
    build_all_satellite_positions,
    build_visibility_matrix,
    gateway_ecef,
    NUM_SATELLITES,
)
from .gateways import (
    GATEWAY_LOCATIONS,
    NUM_GATEWAYS,
    GATEWAY_CAPACITY_GBPS,
)
from .routing import (
    build_correlation_matrix,
    GORA,
    AORA,
)

# Called once per execution
# Builds every possible traffic sending combination
# approx. 2070 pairs & then shuffles randomly
# each execution processes demands in a different order
# earlier demands get first access to link capacity
def _makeAllGatewayPairs(numGateways: int, rng):
    allPairs = [(sourceIdx, destIdx)
                for sourceIdx in range(numGateways)
                for destIdx in range(numGateways)
                if sourceIdx != destIdx]

    # shuffle in place so demands arrive in a random order each execution
    rng.shuffle(allPairs)
    return allPairs

def run_simulation(
    numDemands: int = 200, # how many routing request per execution
    trafficLoadLevels = None, # which percentages to test
                              # default = 13 points spaces 0 to 120%
    randomSeed: int = 42, # identical results every time
    extraSearchRadiusKm: float = 600.0, # AORA's search area size
    satelliteSnapshotTime: float = 0.0, # when to freeze satellite positions
    demandBandwidthGbps: float = 1.0, # bandwidth each request needs
):
    # default to 13 evenly spaced points: 0.0, 0.1, 0.2, ... 1.2
    if trafficLoadLevels is None:
        trafficLoadLevels = np.linspace(0.0, 1.2, 13)

    pythonRng = random.Random(randomSeed)

    # satellite and gateway positions computed once reused across all load levels
    print("Building satellite positions...")
    satellitePositionsEcef = build_all_satellite_positions(satelliteSnapshotTime) # (720, 3)
    gatewayPositionsEcef = np.array(
        [gateway_ecef(lat, lon) for lat, lon in GATEWAY_LOCATIONS]
    )

    print("Computing visibility matrix...")
    visibilityMatrix, slantRangesKm = build_visibility_matrix( # (46, 720) bool and float
        gatewayPositionsEcef, satellitePositionsEcef
    )
    print(f"  Visible gateway-satellite links: {visibilityMatrix.sum()}")

    numGateways = NUM_GATEWAYS # 46
    numSatellites = NUM_SATELLITES # 720
    fullMatrixSize = numGateways + numSatellites # 766 nodes total
    totalNetworkCapacityGbps = numGateways * GATEWAY_CAPACITY_GBPS # 42,320 Gbps
    numLoadLevels = len(trafficLoadLevels)

    # one entry per load level filled in during the main loop below
    goraMatrixSize = np.full(numLoadLevels, fullMatrixSize, dtype=int)
    aoraAvgSubmatrixSize = np.zeros(numLoadLevels)
    goraAvgCompTimeMs = np.zeros(numLoadLevels)
    aoraAvgCompTimeMs = np.zeros(numLoadLevels)
    goraRoutingCapacity = np.zeros(numLoadLevels)
    aoraRoutingCapacity = np.zeros(numLoadLevels)
    goraAvgDelayMs = np.zeros(numLoadLevels)
    aoraAvgDelayMs = np.zeros(numLoadLevels)
    aoraFallbackRate = np.zeros(numLoadLevels)
    aoraSubmatrixSizesPerLevel = [[] for _ in range(numLoadLevels)]

    # Generate a single fixed demand set used at every load level.
    # Keeping the pairs constant isolates load (bandwidth) as the only variable,
    # which lets GORA vs AORA routing quality differences accumulate cleanly.
    fixedGatewayPairs = _makeAllGatewayPairs(numGateways, pythonRng)[:numDemands]

    for loadLevelIdx, currentLoadLevel in enumerate(trafficLoadLevels):
        print(f"\nLoad level {currentLoadLevel:.2f} ...")

        # Always route numDemands demands for statistical validity.
        # Scale per-demand bandwidth so total offered traffic = loadLevel × network capacity.
        # At 0% load, use minimum bandwidth (tests pure topological connectivity).
        if currentLoadLevel > 0.0:
            scaledBandwidthGbps = currentLoadLevel * totalNetworkCapacityGbps / numDemands
        else:
            scaledBandwidthGbps = demandBandwidthGbps  # 1 Gbps default

        demandsThisLevel = fixedGatewayPairs

        # GORA PASS link loads start at zero and grow as demands are accepted
        goraLinkLoadsGbps = np.zeros((numGateways, numSatellites), dtype=float) # (46, 720)

        # full 766x766 graph edge weights are current link loads (inf if no LOS)
        goraCorrelationMatrix = build_correlation_matrix(
            GATEWAY_LOCATIONS, satellitePositionsEcef,
            visibilityMatrix, slantRangesKm,
            goraLinkLoadsGbps
        )

        goraRouter = GORA(
            goraCorrelationMatrix, gatewayPositionsEcef,
            satellitePositionsEcef, numGateways
        )

        goraComputationTimesMs = []
        goraAcceptedDelaysMs = []
        goraNumRouted = 0

        for sourceGwIdx, destGwIdx in demandsThisLevel:
            # dijkstra on the full 766-node matrix
            routedPath, _, _, propagationDelayMs, computeTimeMs = \
                goraRouter.route(sourceGwIdx, destGwIdx)

            goraComputationTimesMs.append(computeTimeMs)

            if routedPath:
                # check every hop, reject demand if any link would exceed 920 Gbps
                demandIsFeasible = True
                pendingHopUpdates = []

                for nodeA, nodeB in zip(routedPath[:-1], routedPath[1:]):
                    # gateways are 0-45; satellites are 46-765 in matrix but 0-719 in link table
                    if nodeA < numGateways and nodeB >= numGateways:
                        gatewayIdx = nodeA
                        satelliteIdx = nodeB - numGateways
                    elif nodeB < numGateways and nodeA >= numGateways:
                        gatewayIdx = nodeB
                        satelliteIdx = nodeA - numGateways
                    else:
                        continue # gw-to-gw hop, no direct links, skip

                    projectedLoad = goraLinkLoadsGbps[gatewayIdx, satelliteIdx] + scaledBandwidthGbps
                    if projectedLoad > GATEWAY_CAPACITY_GBPS:
                        demandIsFeasible = False
                        break # no point checking the rest

                    pendingHopUpdates.append((gatewayIdx, satelliteIdx, projectedLoad)) # stage update

                if demandIsFeasible:
                    goraNumRouted += 1
                    goraAcceptedDelaysMs.append(propagationDelayMs)

                    # apply staged updates and rebuild so future demands see new congestion
                    for gatewayIdx, satelliteIdx, newLoad in pendingHopUpdates:
                        goraLinkLoadsGbps[gatewayIdx, satelliteIdx] = newLoad

                    goraCorrelationMatrix = build_correlation_matrix(
                        GATEWAY_LOCATIONS, satellitePositionsEcef,
                        visibilityMatrix, slantRangesKm,
                        goraLinkLoadsGbps
                    )
                    goraRouter = GORA(
                        goraCorrelationMatrix, gatewayPositionsEcef,
                        satellitePositionsEcef, numGateways
                    )

        # summarise GORA results for this load level
        goraAvgCompTimeMs[loadLevelIdx] = np.mean(goraComputationTimesMs) if goraComputationTimesMs else 0.0
        goraRoutingCapacity[loadLevelIdx] = goraNumRouted / len(demandsThisLevel) if demandsThisLevel else 0.0
        goraAvgDelayMs[loadLevelIdx] = np.mean(goraAcceptedDelaysMs) if goraAcceptedDelaysMs else 0.0

        # AORA PASS, same demands, fresh link-load state for a fair comparison
        aoraLinkLoadsGbps = np.zeros((numGateways, numSatellites), dtype=float)

        aoraCorrelationMatrix = build_correlation_matrix(
            GATEWAY_LOCATIONS, satellitePositionsEcef,
            visibilityMatrix, slantRangesKm,
            aoraLinkLoadsGbps
        )

        aoraRouter = AORA(
            aoraCorrelationMatrix, GATEWAY_LOCATIONS,
            gatewayPositionsEcef, satellitePositionsEcef,
            visibilityMatrix, numGateways, extraSearchRadiusKm
        )

        aoraComputationTimesMs = []
        aoraAcceptedDelaysMs = []
        aoraNumRouted = 0
        aoraSubmatrixSizes = []
        aoraNumFallbacks = 0

        for sourceGwIdx, destGwIdx in demandsThisLevel:
            # also returns subMatrixSize and usedFallback unlike GORA
            routedPath, _, subMatrixSize, propagationDelayMs, computeTimeMs, usedFallback = \
                aoraRouter.route(sourceGwIdx, destGwIdx)

            aoraComputationTimesMs.append(computeTimeMs)
            aoraSubmatrixSizes.append(subMatrixSize)

            if usedFallback:
                aoraNumFallbacks += 1 # sub-matrix had no path, used full 766-node matrix

            if routedPath and not usedFallback:
                # Count only sub-matrix successes toward AORA routing capacity.
                # Fallback demands are excluded so we measure AORA's standalone quality,
                # matching how the paper separates GORA from AORA performance.
                demandIsFeasible = True
                pendingHopUpdates = []

                for nodeA, nodeB in zip(routedPath[:-1], routedPath[1:]):
                    if nodeA < numGateways and nodeB >= numGateways:
                        gatewayIdx = nodeA
                        satelliteIdx = nodeB - numGateways
                    elif nodeB < numGateways and nodeA >= numGateways:
                        gatewayIdx = nodeB
                        satelliteIdx = nodeA - numGateways
                    else:
                        continue

                    projectedLoad = aoraLinkLoadsGbps[gatewayIdx, satelliteIdx] + scaledBandwidthGbps
                    if projectedLoad > GATEWAY_CAPACITY_GBPS:
                        demandIsFeasible = False
                        break

                    pendingHopUpdates.append((gatewayIdx, satelliteIdx, projectedLoad))

                if demandIsFeasible:
                    aoraNumRouted += 1
                    aoraAcceptedDelaysMs.append(propagationDelayMs)

                    for gatewayIdx, satelliteIdx, newLoad in pendingHopUpdates:
                        aoraLinkLoadsGbps[gatewayIdx, satelliteIdx] = newLoad

                    # Rebuild correlation matrix with updated AORA link loads
                    aoraCorrelationMatrix = build_correlation_matrix(
                        GATEWAY_LOCATIONS, satellitePositionsEcef,
                        visibilityMatrix, slantRangesKm,
                        aoraLinkLoadsGbps
                    )
                    aoraRouter = AORA(
                        aoraCorrelationMatrix, GATEWAY_LOCATIONS,
                        gatewayPositionsEcef, satellitePositionsEcef,
                        visibilityMatrix, numGateways, extraSearchRadiusKm
                    )

        # summarise AORA results for this load level
        aoraAvgCompTimeMs[loadLevelIdx] = np.mean(aoraComputationTimesMs) if aoraComputationTimesMs else 0.0
        aoraRoutingCapacity[loadLevelIdx] = aoraNumRouted / len(demandsThisLevel) if demandsThisLevel else 0.0
        aoraAvgDelayMs[loadLevelIdx] = np.mean(aoraAcceptedDelaysMs) if aoraAcceptedDelaysMs else 0.0
        aoraFallbackRate[loadLevelIdx] = aoraNumFallbacks / len(demandsThisLevel) if demandsThisLevel else 0.0
        aoraAvgSubmatrixSize[loadLevelIdx] = np.mean(aoraSubmatrixSizes) if aoraSubmatrixSizes else 0.0

        aoraSubmatrixSizesPerLevel[loadLevelIdx] = aoraSubmatrixSizes # raw sizes for histogram

        print(f"  GORA routed {goraNumRouted}/{len(demandsThisLevel)}, "
              f"avg time {goraAvgCompTimeMs[loadLevelIdx]:.2f} ms, "
              f"avg delay {goraAvgDelayMs[loadLevelIdx]:.1f} ms")
        print(f"  AORA routed {aoraNumRouted}/{len(demandsThisLevel)} (sub-matrix only), "
              f"avg sub-size {aoraAvgSubmatrixSize[loadLevelIdx]:.0f}, "
              f"fallbacks {aoraNumFallbacks}, "
              f"avg time {aoraAvgCompTimeMs[loadLevelIdx]:.2f} ms, "
              f"avg delay {aoraAvgDelayMs[loadLevelIdx]:.1f} ms")

    # pack results and return to main.py for plotting
    resultsDict = {
        'load_levels': trafficLoadLevels,
        'gora_matrix_size': goraMatrixSize,
        'aora_avg_submatrix_size': aoraAvgSubmatrixSize,
        'gora_avg_comp_time_ms': goraAvgCompTimeMs,
        'aora_avg_comp_time_ms': aoraAvgCompTimeMs,
        'gora_routing_capacity': goraRoutingCapacity,
        'aora_routing_capacity': aoraRoutingCapacity,
        'gora_avg_delay_ms': goraAvgDelayMs,
        'aora_avg_delay_ms': aoraAvgDelayMs,
        'aora_fallback_rate': aoraFallbackRate,
        'aora_submatrix_sizes': aoraSubmatrixSizesPerLevel,
        'full_matrix_size': fullMatrixSize,
        'num_gateways': numGateways,
        'num_satellites': numSatellites,
        'extra_radius_km': extraSearchRadiusKm,
        'num_demands': numDemands,
    }
    return resultsDict
