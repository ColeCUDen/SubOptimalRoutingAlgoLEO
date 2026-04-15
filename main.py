import argparse
import numpy as np

from src.simulation import run_simulation
from src.plots import plot_all


# parses command line arguments
def parseArgs():
    p = argparse.ArgumentParser(description="LEO satellite routing simulation")
    p.add_argument("--demands", type=int, default=200, help="routing requests per load level (default 200)")
    p.add_argument("--radius", type=float, default=600.0, help="AORA search radius in km (default 600)")
    p.add_argument("--seed", type=int, default=42, help="random seed (default 42)")
    p.add_argument("--time", type=float, default=0.0, help="satellite snapshot time in seconds (default 0)")
    p.add_argument("--output", type=str, default="figures", help="output directory for figures (default figures)")
    p.add_argument("--loads", type=int, default=13, help="number of load levels in [0%%, 120%%] (default 13)")
    return p.parse_args()


def main():
    args = parseArgs()
    loadLevels = np.linspace(0.0, 1.2, args.loads)

    print(f"demands: {args.demands}, radius: {args.radius} km, seed: {args.seed}, load levels: {args.loads}")

    results = run_simulation(
        numDemands=args.demands,
        trafficLoadLevels=loadLevels,
        randomSeed=args.seed,
        extraSearchRadiusKm=args.radius,
        satelliteSnapshotTime=args.time,
    )

    plot_all(results, outputDir=args.output)

    # summary at midpoint (~60% load)
    mid = len(loadLevels) // 2
    print(f"\nAt ~{loadLevels[mid]*100:.0f}% load:")
    print(f"  GORA matrix size: {results['gora_matrix_size'][mid]}")
    print(f"  AORA avg sub-size: {results['aora_avg_submatrix_size'][mid]:.0f}")
    print(f"  GORA comp time: {results['gora_avg_comp_time_ms'][mid]:.2f} ms")
    print(f"  AORA comp time: {results['aora_avg_comp_time_ms'][mid]:.2f} ms")
    print(f"  GORA routing capacity: {results['gora_routing_capacity'][mid]*100:.1f}%")
    print(f"  AORA routing capacity: {results['aora_routing_capacity'][mid]*100:.1f}%")
    print(f"  GORA avg delay: {results['gora_avg_delay_ms'][mid]:.1f} ms")
    print(f"  AORA avg delay: {results['aora_avg_delay_ms'][mid]:.1f} ms")
    print(f"  AORA fallback rate: {results['aora_fallback_rate'][mid]*100:.1f}%")


if __name__ == "__main__":
    main()
