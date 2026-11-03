"""The three MayaBRIC pillars and the composite (manuscript §3.3).

A  Salt-Source Accessibility       A = exp(-lambda * cost_hr)
D  Regional Administrative Hub     D = 1 - exp(-n_hub25km / k),
   Density                             n_hub25km = sum(1/BWRANK) of ranked EAAMS within 25 km
G  Institutional connectivity      G = 1 - exp(-G_input / a),
                                       G_input = max over the 30 nearest ranked EAAMS of
                                       (1/BWRANK) / max(distance_km, 0.1)
MayaBRIC = (A + D + G) / 3, banded at 0.25 / 0.50 / 0.75.

EAAMS is the Electronic Atlas of Ancient Maya Sites (Witschey & Brown). Only sites
with a numeric BWRANK enter D and G; ranks such as "c" and "sm" are excluded.
"""
import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.spatial import cKDTree

from .constants import (A_G, BAND_CUTS, D_RING_M, G_FLOOR_KM, G_K_NEAREST, K_D,
                        LAMBDA_HR, UTM16N, WGS84)

_to_utm = Transformer.from_crs(WGS84, UTM16N, always_xy=True)


def to_utm(lonlat):
    return np.column_stack(_to_utm.transform(lonlat[:, 0], lonlat[:, 1]))


class RankedSites:
    """Ranked EAAMS sites, projected to UTM 16N, with 1/BWRANK weights."""

    def __init__(self, lonlat, bwrank, names=None):
        rank = pd.to_numeric(pd.Series(bwrank), errors="coerce").values
        keep = ~np.isnan(rank)
        self.xy = to_utm(np.asarray(lonlat, dtype=float)[keep])
        self.rank = rank[keep].astype(int)
        self.weight = 1.0 / self.rank
        self.names = None if names is None else np.asarray(names)[keep]
        self.tree = cKDTree(self.xy)


def n_hub25km(points_xy, sites):
    """Sum of 1/BWRANK over ranked sites within 25 km of each point."""
    return np.array([sites.weight[i].sum() for i in sites.tree.query_ball_point(points_xy, D_RING_M)])


def g_input(points_xy, sites):
    """Strongest single gravity tie: max of (1/BWRANK) / max(km, 0.1) over the 30 nearest.

    Returns (G_input, index into the ranked sites of the winning hub, its true distance in km).
    """
    k = min(G_K_NEAREST, len(sites.rank))
    d_m, idx = sites.tree.query(points_xy, k=k)
    d_km = d_m / 1000.0
    grav = sites.weight[idx] / np.maximum(d_km, G_FLOOR_KM)
    best = grav.argmax(axis=1)
    rows = np.arange(len(points_xy))
    return grav[rows, best], idx[rows, best], d_km[rows, best]


def a_score(cost_hr, lam=LAMBDA_HR):
    return np.exp(-lam * np.asarray(cost_hr, dtype=float))


def d_score(n_hub):
    return 1.0 - np.exp(-np.asarray(n_hub, dtype=float) / K_D)


def g_score(g_in):
    return 1.0 - np.exp(-np.asarray(g_in, dtype=float) / A_G)


def composite(a, d, g):
    return (np.asarray(a) + np.asarray(d) + np.asarray(g)) / 3.0


def band(score):
    """Fixed bands: 0 below 0.25, 1 from 0.25, 2 from 0.50, 3 from 0.75."""
    return np.searchsorted(BAND_CUTS, np.asarray(score), side="right").astype(np.int32)


def median_anchored_lambda(cost_hr):
    """lambda = 1 / median cost, so that a structure at the median cost scores 1/e.

    The published catalog uses this value (0.0424784, reported as 0.04248) for
    every detection; query sites use the reported 0.04248.
    """
    return 1.0 / float(np.median(cost_hr))
