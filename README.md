# LEO Satellite Routing — Liu & Zhu (2018) Reproduction

Python reproduction of **"A Suboptimal Routing Algorithm for Massive LEO Satellite Networks"** (Liu & Zhu, ISNCC 2018).

Implements GORA and AORA for a 720-satellite, 46-gateway LEO constellation and reproduces Figures 4–7 from the paper.

---

## System overview

| Parameter | Value |
|---|---|
| Satellites | 720 (18 planes × 40 per plane) |
| Altitude | 1,200 km |
| Inclination | 90° (polar) |
| Walker phasing | 4.5° between adjacent planes |
| Ground gateways | 46 globally distributed |
| Gateway capacity | 920 Gbps each (42,320 Gbps total) |
| Correlation matrix | 766 × 766 (46 gateways + 720 satellites) |
| Visibility threshold | Elevation angle ≥ 10° |

---

## Repository structure

```
SubOptimalRoutingAlgoLEO/
├── main.py                 # entry point
├── requirements.txt
├── evaluation_draft.md     # methodology, metrics, and preliminary results
├── src/
│   ├── __init__.py
│   ├── constellation.py    # orbital mechanics, satellite positions, visibility
│   ├── gateways.py         # 46 gateway coordinates (lat/lon)
│   ├── routing.py          # GORA, AORA, Dijkstra wrapper, helper functions
│   ├── simulation.py       # main simulation loop
│   └── plots.py            # figure generation (Figs 4–7 + bonus)
├── figures/                # output figures (created at runtime)
└── README.md
```

---

## Dependencies

- Python ≥ 3.8
- numpy ≥ 1.21
- scipy ≥ 1.7
- matplotlib ≥ 3.4

---

## Setup

```bash
# 1. enter the project directory
cd SubOptimalRoutingAlgoLEO

# 2. (optional) create a virtual environment
python -m venv venv
# Windows:
venv\Scripts\activate
# Linux/macOS:
source venv/bin/activate

# 3. install dependencies
pip install -r requirements.txt
```

---

## Running the simulation

### Default run

```bash
python main.py
```

Expected runtime: **2–4 minutes** on a typical laptop.
For a quick preview: `python main.py --demands 50 --loads 7` (~30 seconds).

Figures are saved to `figures/`:
- `fig4_matrix_size.png`
- `fig5_comp_time.png`
- `fig6_routing_capacity.png`
- `fig7_delay.png`
- `bonus_submatrix_distribution.png`

### Command-line options

| Flag | Default | Description |
|---|---|---|
| `--demands N` | 200 | routing requests per traffic-load level |
| `--radius R` | 600.0 | AORA extra search radius (km) |
| `--seed S` | 42 | random seed |
| `--time T` | 0.0 | satellite snapshot time (seconds) |
| `--loads L` | 13 | number of load levels sampled in [0%, 120%] |
| `--output DIR` | figures | output directory for figures |

**Higher fidelity run (~15–30 min):**
```bash
python main.py --demands 500 --loads 25
```

---

## Reproducing paper figures

| Output file | Paper figure | Description |
|---|---|---|
| `fig4_matrix_size.png` | Fig. 4 | correlation matrix node count vs load |
| `fig5_comp_time.png` | Fig. 5 | average Dijkstra computation time vs load |
| `fig6_routing_capacity.png` | Fig. 6 | fraction of demands routed vs load |
| `fig7_delay.png` | Fig. 7 | average end-to-end propagation delay vs load |

Expected qualitative results:
- **Fig 4/5**: AORA sub-matrix is significantly smaller than GORA's full 766-node matrix, with proportionally lower computation time.
- **Fig 6**: Both algorithms route ~100% of demands at low load; GORA maintains higher capacity above ~80% load.
- **Fig 7**: AORA delay is within a few ms of GORA across all load levels, confirming approximate optimality.

---

## Algorithm notes

### Orbital model
Analytical circular orbit with Walker-Delta phasing — no TLE data needed. RAAN is evenly distributed across 18 planes; mean anomaly includes the 4.5° phasing offset.

### Visibility
Elevation angle ≥ 10° computed from ECEF positions. Links below the threshold are set to `np.inf` in the correlation matrix.

### Dijkstra
`scipy.sparse.csgraph.dijkstra` with `directed=False`. scipy treats 0 as no edge, so all valid edges receive a +1e-9 offset before solving.

### AORA sub-matrix construction
1. Compute spherical midpoint *k* of source and destination gateways.
2. Set radius *R* = haversine(*k*, source) + `extra_radius_km`.
3. Include all gateways within *R* of *k* (source/destination always included).
4. Include all satellites visible to any selected gateway.
5. Run Dijkstra on the sub-matrix.
6. Fall back to the full 766×766 matrix if no path is found.

### Traffic model
Demands are sampled uniformly from all 46 × 45 ordered gateway pairs. Each demand requests 1 Gbps. A demand is accepted only if no link along the path would exceed 920 Gbps.
