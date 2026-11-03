"""Score every detection: D and G inputs from EAAMS, A from least-cost hours, MayaBRIC and band.

  python scripts/score_catalog.py --points detections.csv --costs costs.csv \
      --eaams eaams_points.shp --out scored.csv

--costs is the output of cost_to_salt.py (same row order as --points). By default
lambda is anchored to the median cost of the scored set, as in the published
catalog; pass --lambda 0.04248 to use the published value on a different set.
Stored precision follows the catalog: cost_hr and n_hub25km 4 dp, A 5 dp,
D and G 4 dp, MayaBRIC from the rounded pillars.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mayabric import pillars as P  # noqa: E402
from mayabric.io import read_eaams, read_points  # noqa: E402


def score(pts_lonlat, cost_hr, sites, lam=None):
    xy = P.to_utm(pts_lonlat)
    if lam is None:
        lam = P.median_anchored_lambda(cost_hr)
    n25 = np.round(P.n_hub25km(xy, sites), 4)   # D is computed from the stored 4 dp input
    g_in, hub, hub_km = P.g_input(xy, sites)
    A = np.round(P.a_score(cost_hr, lam), 5)
    D = np.round(P.d_score(n25), 4)
    G = np.round(P.g_score(g_in), 4)
    mb = P.composite(A, D, G)
    out = pd.DataFrame({
        "cost_hr": np.round(cost_hr, 4), "n_hub25km": n25, "G_input": g_in,
        "A_score": A, "D_score": D, "G_score": G, "MayaBRIC": mb, "MB_qband": P.band(mb),
        "d_hub_km": np.round(hub_km, 4), "hub_rank": sites.rank[hub],
    })
    if sites.names is not None:
        out["hub_name"] = sites.names[hub]
    return out, lam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", required=True)
    ap.add_argument("--lat", default="cen_lat")
    ap.add_argument("--lon", default="cen_lon")
    ap.add_argument("--costs", required=True, help="CSV from cost_to_salt.py (cost_hr, best_src)")
    ap.add_argument("--eaams", required=True, help="EAAMS points with BWRANK and BWNAME")
    ap.add_argument("--lambda", dest="lam", type=float, default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tab, pts = read_points(a.points, a.lat, a.lon)
    c = pd.read_csv(a.costs)
    assert len(c) == len(tab), "costs and points differ in length"
    assert np.allclose(c[["lon", "lat"]].to_numpy(float), pts, atol=1e-6), "costs not aligned to points"
    _, sites, _ = read_eaams(a.eaams)
    out, lam = score(pts, c.cost_hr.to_numpy(float), sites, a.lam)
    out.insert(1, "best_src", c.best_src.values)
    out.insert(0, "cen_lon", pts[:, 0])
    out.insert(0, "cen_lat", pts[:, 1])
    out.to_csv(a.out, index=False)
    pct = np.round(np.bincount(out.MB_qband, minlength=4) / len(out) * 100, 2)
    print(f"lambda = {lam:.7f} hr^-1; MayaBRIC median {out.MayaBRIC.median():.4f}; bands % {pct.tolist()}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
