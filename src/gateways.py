"""
gateways.py
46 ground gateway locations (lat, lon in degrees) distributed globally.
Each gateway has a capacity of 920 Gbps (total system = 42,320 Gbps).

Source: Liu & Zhu (2018), Table / system description.
"""

# (name, latitude_deg, longitude_deg)
GATEWAY_DATA = [
    # ── North America (10) ────────────────────────────────────────────────────
    ("New York",       40.71, -74.01),
    ("Los Angeles",    34.05, -118.24),
    ("Chicago",        41.88,  -87.63),
    ("Houston",        29.76,  -95.37),
    ("Seattle",        47.61, -122.33),
    ("Miami",          25.77,  -80.19),
    ("Toronto",        43.65,  -79.38),
    ("Montreal",       45.50,  -73.57),
    ("Vancouver",      49.25, -123.12),
    ("Mexico City",    19.43,  -99.13),
    # ── South America (5) ─────────────────────────────────────────────────────
    ("Sao Paulo",     -23.55,  -46.63),
    ("Buenos Aires",  -34.60,  -58.38),
    ("Lima",          -12.05,  -77.04),
    ("Bogota",          4.71,  -74.07),
    ("Santiago",      -33.45,  -70.67),
    # ── Europe (8) ────────────────────────────────────────────────────────────
    ("London",         51.51,   -0.13),
    ("Paris",          48.86,    2.35),
    ("Berlin",         52.52,   13.40),
    ("Madrid",         40.42,   -3.70),
    ("Rome",           41.90,   12.50),
    ("Stockholm",      59.33,   18.07),
    ("Amsterdam",      52.37,    4.90),
    ("Prague",         50.08,   14.44),
    # ── Africa (6) ────────────────────────────────────────────────────────────
    ("Cairo",          30.04,   31.24),
    ("Johannesburg",  -26.20,   28.04),
    ("Lagos",           6.46,    3.39),
    ("Nairobi",        -1.29,   36.82),
    ("Tunis",          36.82,   10.17),
    ("Khartoum",       15.55,   32.53),
    # ── Middle East (3) ───────────────────────────────────────────────────────
    ("Abu Dhabi",      24.47,   54.37),
    ("Baghdad",        33.34,   44.40),
    ("Tehran",         35.69,   51.39),
    # ── Asia (10) ─────────────────────────────────────────────────────────────
    ("Moscow",         55.75,   37.62),
    ("Beijing",        39.91,  116.39),
    ("Shanghai",       31.23,  121.47),
    ("Hong Kong",      22.32,  114.17),
    ("Tokyo",          35.68,  139.69),
    ("Seoul",          37.57,  126.98),
    ("New Delhi",      28.61,   77.21),
    ("Singapore",       1.35,  103.82),
    ("Bangkok",        13.75,  100.52),
    ("Hanoi",          21.03,  105.85),
    # ── Oceania (4) ───────────────────────────────────────────────────────────
    ("Sydney",        -33.87,  151.21),
    ("Melbourne",     -37.81,  144.96),
    ("Auckland",      -36.85,  174.76),
    ("Brisbane",      -27.47,  153.03),
]

NUM_GATEWAYS          = len(GATEWAY_DATA)          # 46
GATEWAY_CAPACITY_GBPS = 920.0                       # per gateway
TOTAL_CAPACITY_GBPS   = NUM_GATEWAYS * GATEWAY_CAPACITY_GBPS  # 42,320 Gbps

# Convenience: list of (lat, lon) tuples
GATEWAY_LOCATIONS = [(lat, lon) for _, lat, lon in GATEWAY_DATA]
GATEWAY_NAMES     = [name       for name, _, _  in GATEWAY_DATA]
