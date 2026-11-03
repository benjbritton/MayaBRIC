"""Site-level clusters (DBSCAN, UTM 16N, eps 750 m, min_samples 5), spiral-numbered,
with each cluster's nearest recorded site and a flag for no site within 10 km.

  python scripts/clusters.py --points detections.csv --eaams eaams_points.shp \
      --out-points cluster_membership.csv --out-clusters clusters.csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from pyproj import Transformer

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mayabric.clusters import dbscan_labels, nearest_site_km, spiral_numbering  # noqa: E402
from mayabric.constants import UTM16N, WGS84  # noqa: E402
from mayabric.io import read_eaams, read_points  # noqa: E402
from mayabric.pillars import to_utm  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", required=True)
    ap.add_argument("--lat", default="cen_lat")
    ap.add_argument("--lon", default="cen_lon")
    ap.add_argument("--eaams", required=True, help="EAAMS points (all sites are used for distance)")
    ap.add_argument("--out-points", required=True)
    ap.add_argument("--out-clusters", required=True)
    a = ap.parse_args()

    _, pts = read_points(a.points, a.lat, a.lon)
    xy = to_utm(pts)
    labels = dbscan_labels(xy)
    cid, csz, t = spiral_numbering(labels, xy)
    site_ll, _, e = read_eaams(a.eaams)
    d_pt, _ = nearest_site_km(pts, site_ll)
    pd.DataFrame({"cen_lat": pts[:, 1], "cen_lon": pts[:, 0], "cluster_id": cid, "clust_sz": csz,
                  "d_eaams_km": np.round(d_pt, 4)}).to_csv(a.out_points, index=False)

    lon, lat = Transformer.from_crs(UTM16N, WGS84, always_xy=True).transform(t[:, 4], t[:, 5])
    cl_ll = np.column_stack([lon, lat])
    d_cl, i_cl = nearest_site_km(cl_ll, site_ll)
    names = e["BWNAME"].where(e["BWNAME"].notna(), "unnamed").astype(str).values
    cl = pd.DataFrame({
        "cluster_id": t[:, 1].astype(int), "dbscan_lbl": t[:, 0].astype(int), "ring": t[:, 2].astype(int),
        "n_struct": t[:, 3].astype(int), "lat": np.round(lat, 6), "lon": np.round(lon, 6),
        "near_site": names[i_cl], "near_rank": e["BWRANK"].astype(str).values[i_cl],
        "near_km": np.round(d_cl, 4), "NoEAAMS10k": (d_cl > 10).astype(int),
    }).sort_values("cluster_id")
    cl.to_csv(a.out_clusters, index=False)
    print(f"{len(cl)} clusters, {(cid == 0).sum():,} unclustered detections, "
          f"{cl.NoEAAMS10k.sum()} clusters with no recorded site within 10 km")


if __name__ == "__main__":
    main()
