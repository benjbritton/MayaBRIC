"""Published MayaBRIC constants (manuscript §3.3; Supplementary S2).

Every pillar maps a raw, survey-portable input onto a bounded score with one
constant anchored to the catalog median of that input. Changing a constant
rescales the score without touching the underlying measurement.
"""

# A pillar: Salt-Source Accessibility. A = exp(-LAMBDA_HR * cost_hr)
LAMBDA_HR = 0.04248          # hr^-1; = 1 / median catalog cost (23.54 hr)

# D pillar: Regional Administrative Hub Density. D = 1 - exp(-n_hub25km / K_D)
K_D = 14.43
D_RING_M = 25_000            # hard 25 km catchment, about one day's round trip on foot

# G pillar: institutional connectivity (gravity). G = 1 - exp(-G_input / A_G)
A_G = 0.0941
G_K_NEAREST = 30             # ranked EAAMS sites considered per structure
G_FLOOR_KM = 0.1             # distance floor, so a site at its own hub stays finite

# Composite and fixed score bands: 0 < 0.25 <= 1 < 0.50 <= 2 < 0.75 <= 3
BAND_CUTS = (0.25, 0.50, 0.75)

# Site-level clusters (DBSCAN in UTM 16N)
DBSCAN_EPS_M = 750
DBSCAN_MIN_SAMPLES = 5
SPIRAL_RINGS = 3
MERIDA_LONLAT = (-89.6237, 20.9674)   # spiral numbering starts from the cluster nearest Merida

# Least-cost surface
TOBLER_FLAT_SLOPE = 0.0      # Tobler multiplier is relative to walking on flat ground
BARRIER = 1e5                # friction >= this value is impassable (permanent water = 1e6)
SITE_LCP_PAD_M = 20_000      # padding around a point pair for the windowed site LCP
SITE_N_SOURCES = 5           # nearest salt positions (straight line) tried per query site

# Coordinate reference systems
WGS84 = "EPSG:4326"
UTM16N = "EPSG:32616"        # distances for D, G and clusters
LCC = ("+proj=lcc +lat_1=17.0 +lat_2=21.0 +lat_0=19.0 +lon_0=-89.0 "
       "+x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")  # cost surface
EARTH_RADIUS_KM = 6371.0
