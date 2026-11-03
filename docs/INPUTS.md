# Inputs

| Input | Used for | Where to get it | In this repository? |
|---|---|---|---|
| Salt-source anchors | A pillar | `data/salt_sources.csv`: the twelve positions of manuscript Table 1 | Yes |
| Detections (structure centroids) | everything | The Q2000 Master Shapefile is restricted to protect archaeological sites. It is available from the corresponding author to identified researchers under the terms in Supplementary S3. Any table of point centroids works: pass `--lat` / `--lon` for the column names | No |
| EAAMS site points (`BWRANK`, `BWNAME`) | D and G pillars, cluster distances | Electronic Atlas of Ancient Maya Sites, W. R. T. Witschey and C. T. Brown. Available from its authors under their terms | No |
| Friction surface, 30 m | A pillar | Built from the public layers below with the parameters in `config/friction_config.yaml` | No (builder to be added) |
| FABDEM 30 m | Slope term of the A pillar | FABDEM (Hawker et al. 2022), CC BY-NC-SA 4.0. Reproject to the LCC grid in `config/friction_config.yaml` | No |
| HydroRIVERS v1.0 | Friction: rivers | Lehner & Grill (2013), hydrosheds.org | No |
| Natural Earth 10 m coastline | Friction: coastal canoe | naturalearthdata.com, public domain | No |
| JRC Global Surface Water, seasonality | Friction: permanent water | Pekel et al. (2016), global-surface-water.appspot.com | No |
| Reprocessed G-LiHT terrain products | Detection (upstream of this code) | OpenTopography (DOI in the manuscript's Data Availability) | No |

The friction raster and the reprojected DEM must share one grid; `cost_to_salt.py`
crops both to their common shape, as the published run did.

EAAMS sites with a non-numeric rank (for example `c`, `sm`) are excluded from D and
G but are kept for the nearest-site distance used to flag clusters with no recorded
site within 10 km.
