"""
main.py
Entry point for the LEO satellite routing simulation.

Reproduces Figures 4–7 from:
  Liu & Zhu (2018) "A Suboptimal Routing Algorithm for Massive LEO
  Satellite Networks", ISNCC 2018.

Usage:
    python main.py                      # default parameters
    python main.py --demands 300        # more demands per load level
    python main.py --radius 800         # larger AORA search radius (km)
    python main.py --seed 7             # different random seed
    python main.py --time 0             # satellite snapshot time (seconds)
    python main.py --output figures     # output directory for plots
"""

import argparse
import numpy as np

from src.simulation import run_simulation
from src.plots import plot_all


def parse_args():
    p = argparse.ArgumentParser(
        description="LEO satellite routing simulation — Liu & Zhu (2018)"
    )
    p.add_argument("--demands", type=int,   default=200,
                   help="Number of random gateway-pair demands per load level (default 200)")
    p.add_argument("--radius",  type=float, default=600.0,
                   help="AORA extra search radius in km (default 600)")
    p.add_argument("--seed",    type=int,   default=42,
                   help="Random seed (default 42)")
    p.add_argument("--time",    type=float, default=0.0,
                   help="Satellite snapshot time in seconds (default 0)")
    p.add_argument("--output",  type=str,   default="figures",
                   help="Output directory for figures (default 'figures')")
    p.add_argument("--loads",   type=int,   default=13,
                   help="Number of traffic load levels sampled in [0, 1.2] (default 13)")
    return p.parse_args()


def main():
    args = parse_args()

    load_levels = np.linspace(0.0, 1.2, args.loads)

    print("=" * 60)
    print("LEO Satellite Routing — Liu & Zhu (2018) Reproduction")
    print("=" * 60)
    print(f"  Demands per level : {args.demands}")
    print(f"  AORA radius       : {args.radius} km")
    print(f"  Random seed       : {args.seed}")
    print(f"  Satellite time    : {args.time} s")
    print(f"  Load levels       : {args.loads} points in [0 %, 120 %]")
    print(f"  Output directory  : {args.output}/")
    print()

    results = run_simulation(
        num_demands         = args.demands,
        traffic_load_levels = load_levels,
        random_seed         = args.seed,
        extra_radius_km     = args.radius,
        time_s              = args.time,
    )

    print("\n" + "=" * 60)
    print("Generating figures …")
    print("=" * 60)
    plot_all(results, output_dir=args.output)

    # ── Summary statistics ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    mid = len(load_levels) // 2   # ~60 % load level index

    print(f"At ~{load_levels[mid]*100:.0f}% traffic load:")
    print(f"  GORA matrix size      : {results['gora_matrix_size'][mid]}")
    print(f"  AORA avg sub-size     : {results['aora_avg_submatrix_size'][mid]:.0f}")
    print(f"  GORA avg comp time    : {results['gora_avg_comp_time_ms'][mid]:.2f} ms")
    print(f"  AORA avg comp time    : {results['aora_avg_comp_time_ms'][mid]:.2f} ms")
    print(f"  GORA routing capacity : {results['gora_routing_capacity'][mid]*100:.1f}%")
    print(f"  AORA routing capacity : {results['aora_routing_capacity'][mid]*100:.1f}%")
    print(f"  GORA avg delay        : {results['gora_avg_delay_ms'][mid]:.1f} ms")
    print(f"  AORA avg delay        : {results['aora_avg_delay_ms'][mid]:.1f} ms")
    print(f"  AORA fallback rate    : {results['aora_fallback_rate'][mid]*100:.1f}%")


if __name__ == "__main__":
    main()
