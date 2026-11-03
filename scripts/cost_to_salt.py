"""Least-cost hours from every detection to its cheapest salt source (A-pillar input).

One full-raster solve per unique source cell; about 9 minutes and ~10 GB RAM per
solve on the 30 m Yucatan grid, so about 2 hours for the twelve anchors.

  python scripts/cost_to_salt.py --points detections.csv --friction friction_A_dry.tif \
      --dem fabdem_lcc_30m.tif --out costs.csv
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mayabric.cost import costs_from_sources, load_cost_surface, nearest_source  # noqa: E402
from mayabric.io import SALT_CSV, read_points, read_salt_sources  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", required=True, help="detections: CSV or vector file with centroid columns")
    ap.add_argument("--lat", default="cen_lat")
    ap.add_argument("--lon", default="cen_lon")
    ap.add_argument("--friction", required=True, help="class-friction GeoTIFF (LCC grid, 30 m)")
    ap.add_argument("--dem", required=True, help="DEM GeoTIFF on the same grid as --friction")
    ap.add_argument("--salt", default=SALT_CSV, help="salt-source anchors (default: data/salt_sources.csv)")
    ap.add_argument("--out", required=True, help="output CSV: lat, lon, cost_hr, best_src")
    a = ap.parse_args()

    t0 = time.time()
    _, pts = read_points(a.points, a.lat, a.lon)
    _, salt_ll, names = read_salt_sources(a.salt)
    print(f"{len(pts):,} points, {len(names)} salt sources; loading cost surface...")
    surface, tr, cell = load_cost_surface(a.friction, a.dem)
    print(f"  grid {surface.shape}, cell {cell:g} m  [{time.time() - t0:.0f} s]")
    m, s2u = costs_from_sources(surface, tr, cell, salt_ll, pts)
    cost, src = nearest_source(m, s2u, names)
    out = pd.DataFrame({"lat": pts[:, 1], "lon": pts[:, 0], "cost_hr": cost, "best_src": src})
    out.to_csv(a.out, index=False)
    print(f"wrote {a.out}: median cost {np.median(cost):.2f} hr  [{time.time() - t0:.0f} s]")
    print(out.best_src.value_counts().to_string())


if __name__ == "__main__":
    main()
