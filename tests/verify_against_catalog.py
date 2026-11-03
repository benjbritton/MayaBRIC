"""Check that this code reproduces the published catalog fields.

Needs the restricted inputs (the Master Shapefile CSV and EAAMS), so it runs only
where those are held. Writes nothing.

  python tests/verify_against_catalog.py --catalog Q2000_master_MB6.csv --eaams eaams_points.shp \
      [--costs costs.csv]

Without --costs the A pillar is checked from the catalog's own cost_hr. With
--costs (output of scripts/cost_to_salt.py) the least-cost hours and the
assigned source are checked as well.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mayabric import pillars as P  # noqa: E402
from mayabric.clusters import dbscan_labels, nearest_site_km, spiral_numbering  # noqa: E402
from mayabric.io import read_eaams  # noqa: E402

FAILED = []


def check(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILED.append(name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--eaams", required=True)
    ap.add_argument("--costs", default=None)
    a = ap.parse_args()

    m = pd.read_csv(a.catalog, keep_default_na=False, na_values=[""])
    ll = m[["cen_lon", "cen_lat"]].to_numpy(float)
    site_ll, sites, _ = read_eaams(a.eaams)
    print(f"{len(m):,} detections; {len(site_ll):,} EAAMS sites, {len(sites.rank):,} ranked")
    xy = P.to_utm(ll)

    print("\nD and G inputs")
    n25 = P.n_hub25km(xy, sites)
    d = np.abs(n25 - m.n_hub25km.values).max()
    check("n_hub25km", d <= 5e-5, f"max |diff| {d:.2e} (catalog stores 4 dp)")
    g_in, hub, hub_km = P.g_input(xy, sites)
    r = np.abs(g_in / m.G_input.values - 1).max()
    check("G_input", r <= 1e-9, f"max relative diff {r:.2e}")
    check("hub_rank", (sites.rank[hub] == m.hub_rank.values).all(), "winning hub rank identical")
    d = np.abs(np.round(hub_km, 4) - m.d_hub_km.values).max()
    check("d_hub_km", d <= 1e-4, f"max |diff| {d:.2e}")

    print("\nPillars, composite, band")
    lam = P.median_anchored_lambda(m.cost_hr.values)
    check("lambda", abs(lam - 0.04248) < 5e-6, f"1/median(cost_hr) = {lam:.7f} (published 0.04248)")
    A = np.round(P.a_score(m.cost_hr.values, lam), 5)
    D = np.round(P.d_score(np.round(n25, 4)), 4)   # catalog computes D from n_hub25km at 4 dp
    G = np.round(P.g_score(g_in), 4)
    for nm, x, tol in (("A_score", A, 1.01e-5), ("D_score", D, 0.0), ("G_score", G, 0.0)):
        diff = np.abs(x - m[nm].values)
        check(nm, diff.max() <= tol, f"max |diff| {diff.max():.1e}; identical in "
              f"{(diff < 1e-12).mean() * 100:.3f}% of rows")
    mb = P.composite(A, D, G)
    diff = np.abs(mb - m.MayaBRIC.values)
    check("MayaBRIC", diff.max() <= 1e-4, f"max |diff| {diff.max():.1e}")
    bd = P.band(mb)
    nb = (bd != m.MB_qband.values).sum()
    check("MB_qband", nb == 0, f"{nb} detections in a different band")
    pct = np.round(np.bincount(bd, minlength=4) / len(bd) * 100, 2).tolist()
    check("band shares", pct == [20.69, 30.95, 43.71, 4.65], f"{pct} % (published 20.69/30.95/43.71/4.65)")

    print("\nClusters")
    labels = dbscan_labels(xy)
    cid, csz, t = spiral_numbering(labels, xy)
    check("cluster count", len(t) == 378 and (cid == 0).sum() == 393,
          f"{len(t)} clusters, {(cid == 0).sum()} unclustered (published 378 / 393)")
    check("cluster_id", (cid == m.cluster_id.values).all(), "identical for every detection")
    check("clust_sz", (csz == m.clust_sz.values).all(), "identical for every detection")
    d_pt, _ = nearest_site_km(ll, site_ll)
    d = np.abs(np.round(d_pt, 4) - m.d_eaams_km.values).max()
    check("d_eaams_km", d <= 1e-4, f"max |diff| {d:.2e}")
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(P.UTM16N, P.WGS84, always_xy=True).transform(t[:, 4], t[:, 5])
    d_cl, _ = nearest_site_km(np.column_stack([lon, lat]), site_ll)
    check("clusters > 10 km from a site", (d_cl > 10).sum() == 166, f"{(d_cl > 10).sum()} (published 166)")

    if a.costs:
        print("\nLeast-cost hours")
        c = pd.read_csv(a.costs)
        diff = np.abs(c.cost_hr.values - m.cost_hr.values)
        check("cost_hr", diff.max() <= 1e-3, f"max |diff| {diff.max():.2e} hr; median "
              f"{np.median(c.cost_hr):.2f} hr (published 23.54)")
        same = (c.best_src.astype(str).values == m.best_src.astype(str).values)
        check("best_src", same.all(), f"identical in {same.mean() * 100:.3f}% of rows")

    print("\n" + ("ALL CHECKS PASSED" if not FAILED else f"FAILED: {', '.join(FAILED)}"))
    sys.exit(1 if FAILED else 0)


if __name__ == "__main__":
    main()
