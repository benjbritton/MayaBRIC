"""Site-level clusters of detections and their distance to recorded sites (§3.2, §5.4).

DBSCAN on structure centroids in UTM 16N (eps 750 m, min_samples 5). Clusters are
numbered 1..N by a clockwise, inward spiral of three distance rings around the
cluster-centroid mean, starting from the cluster nearest Merida; 0 means a
detection belongs to no cluster.
"""
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import BallTree

from .constants import (DBSCAN_EPS_M, DBSCAN_MIN_SAMPLES, EARTH_RADIUS_KM, MERIDA_LONLAT,
                        SPIRAL_RINGS)
from .pillars import to_utm


def dbscan_labels(points_xy):
    return DBSCAN(eps=DBSCAN_EPS_M, min_samples=DBSCAN_MIN_SAMPLES,
                  algorithm="ball_tree").fit(points_xy).labels_


def spiral_numbering(labels, points_xy):
    """Renumber DBSCAN labels by the spiral rule.

    Returns (cluster_id per point, cluster size per point, table of
    raw label, spiral number, ring 1..3, size, centroid x, centroid y).
    """
    ids = np.unique(labels[labels >= 0])
    sizes = np.array([(labels == i).sum() for i in ids])
    cent = np.array([points_xy[labels == i].mean(axis=0) for i in ids])
    merida = to_utm(np.array([MERIDA_LONLAT]))[0]
    centre = cent.mean(axis=0)
    bearing = np.degrees(np.arctan2(cent[:, 0] - centre[0], cent[:, 1] - centre[1])) % 360
    b_start = np.degrees(np.arctan2(*(merida - centre))) % 360
    clockwise = (bearing - b_start) % 360
    radius = np.hypot(*(cent - centre).T)
    cuts = np.quantile(radius, np.linspace(0, 1, SPIRAL_RINGS + 1)[1:-1])
    ring = (SPIRAL_RINGS - 1) - np.searchsorted(cuts, radius, side="right")
    order = np.lexsort((clockwise, ring))
    spiral = np.empty(len(ids), int)
    spiral[order] = np.arange(1, len(ids) + 1)
    num = dict(zip(ids, spiral))
    size = dict(zip(ids, sizes))
    cid = np.array([num.get(l, 0) for l in labels], dtype=np.int32)
    csz = np.array([size.get(l, 0) for l in labels], dtype=np.int32)
    table = np.column_stack([ids, spiral, ring + 1, sizes, cent])
    return cid, csz, table


def nearest_site_km(lonlat, site_lonlat):
    """Great-circle distance (km, R = 6371) from each point to the nearest site.

    Returns (distance km, index of the nearest site).
    """
    tree = BallTree(np.radians(site_lonlat[:, ::-1]), metric="haversine")
    d, i = tree.query(np.radians(lonlat[:, ::-1]), k=1)
    return d[:, 0] * EARTH_RADIUS_KM, i[:, 0]
