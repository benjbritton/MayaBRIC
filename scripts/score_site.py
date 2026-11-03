"""MayaBRIC for named query points outside the detection catalog.

  python scripts/score_site.py --friction friction_A_dry.tif --dem fabdem_lcc_30m.tif \
      --eaams eaams_points.shp "Calakmul Str. II|18.1046|-89.8148" "Uxmal|20.36031|-89.77066"

Each point is "Name|lat|lon". Use --validate to first reproduce the published
Calakmul Structure II to Chetumal Bay cost (49.84 hr).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mayabric.cost import pair_cost_hours  # noqa: E402
from mayabric.io import SALT_CSV, read_eaams, read_salt_sources  # noqa: E402
from mayabric.sites import score_site  # noqa: E402

CALAKMUL = (-89.8148, 18.1046)
CHETUMAL = (-88.300, 18.500)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("points", nargs="*", help='"Name|lat|lon"')
    ap.add_argument("--friction", required=True)
    ap.add_argument("--dem", required=True)
    ap.add_argument("--eaams", required=True)
    ap.add_argument("--salt", default=SALT_CSV)
    ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()

    if a.validate:
        v = pair_cost_hours(a.friction, a.dem, CALAKMUL, CHETUMAL)
        print(f"VALIDATION Calakmul Str. II -> Chetumal Bay: {v:.2f} hr (published 49.84)")
    _, salt_ll, names = read_salt_sources(a.salt)
    _, sites, _ = read_eaams(a.eaams)
    for p in a.points:
        name, lat, lon = p.split("|")
        print(f"\n{name}  {float(lat):.5f}N {abs(float(lon)):.5f}W")
        r = score_site((float(lon), float(lat)), salt_ll, names, sites, a.friction, a.dem)
        print(f"  A = {r['A']:.3f}  (cost {r['cost_hr']:.2f} hr to {r['best_src']})")
        print(f"  D = {r['D']:.3f}  (n_hub25km = {r['n_hub25km']:.3f})")
        print(f"  G = {r['G']:.3f}  (G_input = {r['G_input']:.3f}; hub {r['hub_name']} "
              f"rank {r['hub_rank']} at {r['hub_km']:.2f} km)")
        print(f"  MayaBRIC = {r['MayaBRIC']:.3f}  band {r['band']}")


if __name__ == "__main__":
    main()
