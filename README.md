# Suboptimal Routing Algorithm for Massive LEO Satellite Networks Reproduction

Python reproduction of "A Suboptimal Routing Algorithm for Massive LEO Satellite Networks" (Liu & Zhu, 2018).

Implements two routing algorithms (GORA and AORA) on a 720-satellite, 46-gateway low earth orbit constellation and reproduces figures 4-7 from the paper.

---

## What Is This

The paper proposes a way to route internet traffic through a massive LEO satellite network without having direct satellite-to-satellite links. All traffic hops through ground gateways. GORA finds the globally optimal path by searching the full network. AORA approximates it faster by only searching a smaller region around the midpoint between source and destination.

The simulation runs both algorithms across a range of traffic load levels (0% to 120% of network capacity) and measures how they compare on speed, routing success rate, and signal delay.

---

## System Parameters

- 720 satellites across 18 polar orbital planes (40 per plane) at 1,200 km altitude
- 46 ground gateways distributed globally, 920 Gbps capacity each (42,320 Gbps total)
- Network graph is 766 nodes (46 gateways + 720 satellites)
- A gateway-satellite link exists if the satellite is above 10 degrees elevation

---

## Dependencies

- Python 3.8+
- numpy, scipy, matplotlib

---

## Setup

```bash
cd SubOptimalRoutingAlgoLEO

# optional virtual environment
python -m venv venv
venv\Scripts\activate       # Windows
source venv/bin/activate    # Linux/macOS

pip install -r requirements.txt
```

---

## Running

```bash
python main.py
```

Runs in about 2-4 minutes. Figures save to `figures/`.

Quick preview (~30 seconds):
```bash
python main.py --demands 50 --loads 7
```

Higher fidelity (~15-30 min):
```bash
python main.py --demands 500 --loads 25
```

### Flags

- `--demands` — routing requests per load level (default 200)
- `--radius` — AORA search radius in km (default 600)
- `--seed` — random seed (default 42)
- `--time` — satellite snapshot time in seconds (default 0)
- `--loads` — number of load levels between 0% and 120% (default 13)
- `--output` — output directory for figures (default figures)

---

## Output Figures

- `fig4_matrix_size.png` — matrix node count vs traffic load
- `fig5_comp_time.png` — Dijkstra computation time vs traffic load
- `fig6_routing_capacity.png` — % of demands routed vs traffic load
- `fig7_delay.png` — avg propagation delay vs traffic load
- `bonus_submatrix_distribution.png` — histogram of AORA sub-matrix sizes

---

## How It Works

**GORA** runs Dijkstra on the full 766-node network graph every time. Always finds the optimal path but gets slower as the network grows.

**AORA** finds the midpoint between source and destination, draws a search radius around it, builds a smaller sub-graph from only the nearby gateways and their visible satellites, then runs Dijkstra on that. If no path is found it falls back to GORA. Faster on average but occasionally suboptimal.

Both algorithms update link loads after each accepted demand so future routing decisions reflect current congestion.
