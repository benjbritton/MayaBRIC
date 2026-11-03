"""Slope-adjusted least-cost travel time to salt sources (A-pillar input).

The cost surface is class friction (overland 1.00, river 0.75, coastal canoe
0.70, permanent water barrier 1e6) multiplied by a Tobler slope penalty from a
30 m DEM. Accumulated cost is in seconds per traveller and is reported in hours.

Cost is isotropic (the slope term uses the unsigned gradient magnitude), so the
cost from a source to a structure equals the cost from that structure to the
source. That is what lets one full-raster solve per source serve every
structure at once.
"""
import gc

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds
from skimage.graph import MCP_Geometric

from .constants import BARRIER, LCC, SITE_LCP_PAD_M, WGS84


def tobler_speed(slope):
    """Tobler's hiking function, km/h, for slope as rise/run."""
    return 6.0 * np.exp(-3.5 * np.abs(slope + 0.05))


def slope_adjust(friction, dem, cell_m):
    """Multiply class friction by the Tobler slope penalty; barriers stay barriers."""
    gy, gx = np.gradient(dem.astype(np.float32), cell_m, cell_m)
    slope = np.hypot(gx, gy).astype(np.float32)
    del gx, gy
    mult = tobler_speed(0.0) / tobler_speed(slope)
    del slope
    out = np.where(friction >= BARRIER, friction, friction * mult).astype(np.float32)
    return out


def load_cost_surface(friction_tif, dem_tif, window_bounds=None):
    """Read friction and DEM (same grid), crop to the common shape, slope-adjust.

    Returns (cost_surface, affine transform, cell size in metres).
    """
    with rasterio.open(friction_tif) as f:
        if window_bounds is None:
            fr = f.read(1).astype(np.float32)
            tr = f.transform
        else:
            w = from_bounds(*window_bounds, transform=f.transform)
            fr = f.read(1, window=w).astype(np.float32)
            tr = f.window_transform(w)
    with rasterio.open(dem_tif) as d:
        if window_bounds is None:
            dem = d.read(1).astype(np.float32)
        else:
            dem = d.read(1, window=from_bounds(*window_bounds, transform=d.transform)).astype(np.float32)
    h, w_ = min(fr.shape[0], dem.shape[0]), min(fr.shape[1], dem.shape[1])
    cell_m = abs(tr.a)
    surface = slope_adjust(fr[:h, :w_], dem[:h, :w_], cell_m)
    del fr, dem
    gc.collect()
    return surface, tr, cell_m


def rowcol(x, y, transform, shape):
    """Map projected x, y to a (row, col) cell, clamped to the grid."""
    h, w = shape
    col = int((x - transform.c) / transform.a)
    row = int((y - transform.f) / transform.e)
    return max(0, min(h - 1, row)), max(0, min(w - 1, col))


def costs_from_sources(surface, transform, cell_m, src_lonlat, pts_lonlat, log=print):
    """Accumulated cost (hours) from each source to each point.

    Sources that fall in the same cell are solved once. Returns
    (cost matrix [n_unique_sources, n_points], index of each input source's
    unique cell).
    """
    to_lcc = Transformer.from_crs(WGS84, LCC, always_xy=True)
    sx, sy = to_lcc.transform(src_lonlat[:, 0], src_lonlat[:, 1])
    src_rc = [rowcol(x, y, transform, surface.shape) for x, y in zip(sx, sy)]
    unique_rc = list(dict.fromkeys(src_rc))          # first-appearance order
    for (r, c), x, y in zip(src_rc, src_lonlat[:, 0], src_lonlat[:, 1]):
        if surface[r, c] >= BARRIER:
            log(f"  WARNING: source at {y:.3f}, {x:.3f} lies on a barrier cell; "
                f"every route from it crosses barrier friction, so it is never the cheapest source")
    px, py = to_lcc.transform(pts_lonlat[:, 0], pts_lonlat[:, 1])
    rc = np.array([rowcol(x, y, transform, surface.shape) for x, y in zip(px, py)])
    out = np.full((len(unique_rc), len(rc)), np.inf, dtype=np.float32)
    for i, start in enumerate(unique_rc):
        mcp = MCP_Geometric(surface, sampling=(cell_m, cell_m), fully_connected=True)
        cum, _ = mcp.find_costs(starts=[start])
        out[i] = cum[rc[:, 0], rc[:, 1]] / 3600.0
        del mcp, cum
        gc.collect()
        log(f"  source {i + 1}/{len(unique_rc)} solved")
    return out, np.array([unique_rc.index(r) for r in src_rc])


def nearest_source(cost_matrix, src_to_unique, names):
    """Minimum cost per point and the name of the source that achieves it.

    Where several named sources share one cell, the first listed name is used.
    """
    best_u = cost_matrix.argmin(axis=0)
    first_name = {}
    for i, u in enumerate(src_to_unique):
        first_name.setdefault(u, names[i])
    return cost_matrix.min(axis=0), np.array([first_name[u] for u in best_u])


def pair_cost_hours(friction_tif, dem_tif, a_lonlat, b_lonlat, pad_m=SITE_LCP_PAD_M):
    """Least-cost hours between two points, solved in a padded window.

    Used for scoring individual query sites (manuscript §5.3). A route that would
    need a detour wider than the padding is not available in the window.
    """
    to_lcc = Transformer.from_crs(WGS84, LCC, always_xy=True)
    ax, ay = to_lcc.transform(*a_lonlat)
    bx, by = to_lcc.transform(*b_lonlat)
    bounds = (min(ax, bx) - pad_m, min(ay, by) - pad_m, max(ax, bx) + pad_m, max(ay, by) + pad_m)
    surface, tr, cell_m = load_cost_surface(friction_tif, dem_tif, bounds)
    ra = rowcol(ax, ay, tr, surface.shape)
    rb = rowcol(bx, by, tr, surface.shape)
    cum, _ = MCP_Geometric(surface, sampling=(cell_m, cell_m), fully_connected=True).find_costs(starts=[ra])
    return float(cum[rb] / 3600.0)
