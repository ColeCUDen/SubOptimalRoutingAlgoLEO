# LEO Satellite Routing — Liu & Zhu (2018) Reproduction

Python reproduction of **"A Suboptimal Routing Algorithm for Massive LEO Satellite Networks"** (Liu & Zhu, ISNCC 2018).

Implements the Global Optimal Routing Algorithm (GORA) and the Approximate Optimal Routing Algorithm (AORA) for a 720-satellite, 46-gateway LEO constellation, and reproduces Figures 4–7 from the paper.

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
| Correlation matrix | 766 × 766 (M + N = 46 + 720) |
| Visibility threshold | Elevation angle ≥ 10° |

---

## Repository structure

```
SubOptimalRoutingAlgoLEO/
├── main.py                # Entry point
├── requirements.txt
├── src/
│   ├── __init__.py
│   ├── constellation.py   # Orbital mechanics, satellite positions, visibility
│   ├── gateways.py        # 46 gateway coordinates (lat/lon)
│   ├── routing.py         # GORA, AORA, Dijkstra wrapper, helper functions
│   ├── simulation.py      # Main simulation loop
│   └── plots.py           # Figure generation (Figs 4–7 + bonus)
├── figures/               # Output figures (created at runtime)
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
# 1. Clone / enter the project directory
cd SubOptimalRoutingAlgoLEO

# 2. (Optional) create a virtual environment
python -m venv venv
# Windows:
venv\Scripts\activate
# Linux/macOS:
source venv/bin/activate

# 3. Install dependencies
pip install -r requirements.txt
```

---

## Running the simulation

### Default run (reproduces all figures)

```bash
python main.py
```

Figures are saved to `figures/`:
- `fig4_matrix_size.png`
- `fig5_comp_time.png`
- `fig6_routing_capacity.png`
- `fig7_delay.png`
- `bonus_submatrix_distribution.png`

### Command-line options

| Flag | Default | Description |
|---|---|---|
| `--demands N` | 200 | Demands per traffic-load level |
| `--radius R` | 600.0 | AORA extra search radius (km) |
| `--seed S` | 42 | Random seed |
| `--time T` | 0.0 | Satellite snapshot time (seconds) |
| `--loads L` | 13 | Number of load levels in [0 %, 120 %] |
| `--output DIR` | figures | Output directory |

**Example — higher fidelity run:**
```bash
python main.py --demands 500 --loads 25 --radius 700
```

---

## Reproducing main paper figures

The figures directly correspond to the paper as follows:

| Output file | Paper figure | Description |
|---|---|---|
| `fig4_matrix_size.png` | Fig. 4 | Correlation matrix node count vs load |
| `fig5_comp_time.png` | Fig. 5 | Average Dijkstra computation time vs load |
| `fig6_routing_capacity.png` | Fig. 6 | Fraction of demands routed vs load |
| `fig7_delay.png` | Fig. 7 | Average end-to-end propagation delay vs load |

Expected qualitative results:
- **Fig 4/5**: AORA sub-matrix is significantly smaller than GORA's full 766-node matrix, with a proportionally lower computation time.
- **Fig 6**: Both algorithms route ~100% of demands at low load; GORA maintains higher capacity as load increases above ~80%.
- **Fig 7**: AORA delay is comparable to GORA (within a few ms) across all load levels, confirming "approximate optimality".

---

## Algorithm design notes

### Orbital model
Analytical circular orbit with Walker-Delta phasing — no TLE data needed. RAAN is evenly distributed across 18 planes; mean anomaly uses the phasing offset of 4.5°.

### Visibility
Elevation angle ≥ 10° computed from ECEF Cartesian positions. Links that fail the threshold are set to `np.inf` in the correlation matrix.

### Dijkstra implementation
`scipy.sparse.csgraph.dijkstra` is used with `directed=False`. Because scipy treats weight 0 as "no edge", all valid edges receive a +1e-9 epsilon offset before solving; the offset is subtracted from the returned path cost.

### AORA sub-matrix construction
1. Compute spherical midpoint *k* of source and destination gateways.
2. Radius *R* = haversine(*k*, source) + `extra_radius_km`.
3. Include all gateways within *R* of *k* (source/destination always included).
4. Include all satellites with line-of-sight to any selected gateway.
5. Extract the sub-correlation matrix; run Dijkstra.
6. If no path is found, fall back to the full 766×766 matrix (GORA).

### Traffic model
Demands are sampled uniformly at random from all 46 × 45 ordered gateway pairs. Each demand requests 1 Gbps. A demand is accepted only if no link along the chosen path would exceed 920 Gbps total load.
