"""MayaBRIC for an independently specified location (manuscript §5.1, §5.3).

A query point (for example a site's highest-rank EAAMS entry) is scored with the
same three pillars as a catalog detection. Its A-pillar cost is the cheapest of
windowed least-cost solves to the five salt positions nearest in a straight line.

A point scored at its own EAAMS entry always gets G = 1, and in dense areas
D approaches 1, so query-site scores differ mainly through A.
"""
import numpy as np
from pyproj import Transformer

from .constants import LAMBDA_HR, LCC, SITE_N_SOURCES, WGS84
from .cost import pair_cost_hours
from .pillars import a_score, band, composite, d_score, g_input, g_score, n_hub25km, to_utm

UNREACHABLE_HR = 1e4   # one 30 m barrier cell alone adds about 8,000 hr


def score_site(lonlat, salt_lonlat, salt_names, sites, friction_tif, dem_tif,
               n_sources=SITE_N_SOURCES, lam=LAMBDA_HR, log=print):
    """Return a dict with the cost, the three pillars, MayaBRIC and its band."""
    to_lcc = Transformer.from_crs(WGS84, LCC, always_xy=True)
    q = np.array(to_lcc.transform(*lonlat))
    s = np.column_stack(to_lcc.transform(salt_lonlat[:, 0], salt_lonlat[:, 1]))
    straight_km = np.hypot(*(s - q).T) / 1000.0
    tried = []
    for i in np.argsort(straight_km)[:n_sources]:
        c = pair_cost_hours(friction_tif, dem_tif, tuple(lonlat), tuple(salt_lonlat[i]))
        tried.append((c, salt_names[i], straight_km[i]))
        note = "  (route crosses barrier friction)" if c >= UNREACHABLE_HR else ""
        log(f"  LCP to {salt_names[i]:<30s} {straight_km[i]:6.1f} km straight  {c:7.2f} hr{note}")
    cost, src, _ = min(tried)
    xy = to_utm(np.array([lonlat]))
    n25 = n_hub25km(xy, sites)[0]
    g_in, hub, hub_km = g_input(xy, sites)
    A, D, G = float(a_score(cost, lam)), float(d_score(n25)), float(g_score(g_in[0]))
    mb = float(composite(A, D, G))
    return {
        "cost_hr": cost, "best_src": src, "A": A,
        "n_hub25km": n25, "D": D,
        "G_input": float(g_in[0]), "G": G,
        "hub_name": None if sites.names is None else sites.names[hub[0]],
        "hub_rank": int(sites.rank[hub[0]]), "hub_km": float(hub_km[0]),
        "MayaBRIC": mb, "band": int(band(mb)),
        "tried": tried,
    }
