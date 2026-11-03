"""Readers for the three inputs: salt-source positions, EAAMS sites, detections."""
import os

import geopandas as gpd
import numpy as np
import pandas as pd

from .pillars import RankedSites

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
SALT_CSV = os.path.join(DATA_DIR, "salt_sources.csv")


def read_salt_sources(path=SALT_CSV):
    """Salt-source anchors: columns name, lat, lon (the twelve positions of Table 1)."""
    s = pd.read_csv(path)
    return s, s[["lon", "lat"]].to_numpy(float), s["name"].to_numpy()


def read_eaams(path):
    """EAAMS points (any vector format GeoPandas reads) with BWRANK and BWNAME fields.

    Returns (all-site lon/lat array, RankedSites for numeric ranks only, GeoDataFrame).
    """
    e = gpd.read_file(path).to_crs("EPSG:4326")
    lonlat = np.column_stack([e.geometry.x.values, e.geometry.y.values])
    names = e["BWNAME"].astype(str).values if "BWNAME" in e.columns else None
    return lonlat, RankedSites(lonlat, e["BWRANK"].values, names), e


def read_points(path, lat="cen_lat", lon="cen_lon"):
    """Detections as a table with centroid latitude and longitude columns (CSV or vector file)."""
    if path.lower().endswith(".csv"):
        t = pd.read_csv(path, keep_default_na=False, na_values=[""])
    else:
        t = gpd.read_file(path)
    return t, t[[lon, lat]].to_numpy(float)
