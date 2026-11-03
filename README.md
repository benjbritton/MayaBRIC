# MayaBRIC

Reference implementation of MayaBRIC, the structure-level salt-trade resilience
index in Britton et al., *Mapping 57,181 Maya Structures with LiDAR: From AI Object
Detection to Attribute-Enriched Geospatial Inventory* (Applied Computing and
Geosciences, submitted). It computes, for every detected structure or any query
point:

| Pillar | Input | Score |
|---|---|---|
| **A** Salt-Source Accessibility | slope-adjusted least-cost hours to the cheapest of twelve salt sources | `exp(-0.04248 * cost_hr)` |
| **D** Regional Administrative Hub Density | `n_hub25km` = Σ 1/BWRANK of ranked EAAMS sites within 25 km | `1 - exp(-n_hub25km / 14.43)` |
| **G** Institutional connectivity | `G_input` = max over the 30 nearest ranked EAAMS sites of (1/BWRANK) / max(km, 0.1) | `1 - exp(-G_input / 0.0941)` |

`MayaBRIC = (A + D + G) / 3`, banded at 0.25 / 0.50 / 0.75. Each constant is
anchored to the catalog median of its input (manuscript §3.3; Supplementary S2);
for A, λ = 1/median cost, so a structure at the median cost scores 1/e. The raw
inputs are kept alongside the scores, so the scores can be recalibrated without
recomputing them.

It also produces the site-level clusters of §3.2 (DBSCAN, UTM 16N, eps 750 m,
min_samples 5, spiral-numbered) and each cluster's distance to the nearest recorded
site.

## Reproduces the published catalog

`tests/verify_against_catalog.py` recomputes every MayaBRIC field from the
detection centroids and EAAMS and compares them with the Master Shapefile (MB6,
57,181 detections). Result on 2026-11-02:

| Field | Agreement |
|---|---|
| `n_hub25km`, `G_input`, `hub_rank`, `d_hub_km` | identical to stored precision |
| `D_score`, `G_score` | identical in every row |
| `A_score` | within 1×10⁻⁵ (the catalog's λ comes from unrounded costs; it stores 4 dp) |
| `MayaBRIC` | within 4×10⁻⁵; **band identical for every detection** (20.69 / 30.95 / 43.71 / 4.65 %) |
| `cluster_id`, `clust_sz` | identical for every detection (378 clusters, 393 unclustered) |
| `d_eaams_km`; clusters > 10 km from a recorded site | identical; 166 |
| `cost_hr` from `cost_to_salt.py` | within 4×10⁻⁴ hr; median 23.54 hr (published 23.54) |
| `best_src` from `cost_to_salt.py` | identical for every detection |

The test needs the restricted inputs, so it runs only where they are held.

## Install

```
conda env create -f environment.yml
conda activate mayabric
```

## Use

```
# 1. Least-cost hours to salt (about 2 h for twelve sources on the 30 m Yucatan grid)
python scripts/cost_to_salt.py --points detections.csv \
    --friction friction_A_dry.tif --dem fabdem_lcc_30m.tif --out costs.csv

# 2. Pillars, composite and band for every detection
python scripts/score_catalog.py --points detections.csv --costs costs.csv \
    --eaams eaams_points.shp --out scored.csv

# 3. Site-level clusters
python scripts/clusters.py --points detections.csv --eaams eaams_points.shp \
    --out-points membership.csv --out-clusters clusters.csv

# Query sites outside the catalog (manuscript §5.3), checking Calakmul first
python scripts/score_site.py --validate --friction friction_A_dry.tif \
    --dem fabdem_lcc_30m.tif --eaams eaams_points.shp "Calakmul Str. II|18.1046|-89.8148"
```

`detections.csv` needs centroid columns (`cen_lat`, `cen_lon` by default). See
[docs/INPUTS.md](docs/INPUTS.md) for every input and where to obtain it.

## Notes on the method

- **Cost is isotropic.** Slope enters as the unsigned gradient magnitude, so travel
  costs the same in both directions, and one full-raster solve per source serves all
  structures. Upstream and downstream river travel therefore cost the same.
- **Sources sharing a 30 m cell are solved once.** The twelve anchors fall in twelve
  cells.
- **Query sites** (`score_site.py`) use the reported λ = 0.04248 and a solve windowed
  20 km around each site–source pair, to the five sources nearest in a straight line.
  A point scored at its own EAAMS entry always gets G = 1.
- **EAAMS coverage ends at 17.50° N.** South of it, D and G are not valid.
- **Coastal anchors.** Salt sources are anchored at their documented coastal positions.
  Where an anchor falls on water mapped as impassable in the dry-season friction surface
  (Paynes Creek, Frenchman's Cay), `cost_to_salt.py` reports it. A sensitivity test that
  re-solved both anchors with the surrounding water (10 km and 20 km radii) set to coastal
  canoe friction (0.70) changed the cheapest source for none of the 57,181 detections and
  left Calakmul's two lowest-cost sources (Chetumal Bay, Campeche Bay) unchanged.

## Not yet included

- The friction-surface builder. Its parameters are recorded in
  `config/friction_config.yaml`.
- The wet-season variants and sensitivity analyses of Supplementary S2.8.

## Licence and citation

MIT. Please cite the article and this software; see `CITATION.cff`. The
preprocessing that produced the terrain products is released separately as
[Replicalm](https://github.com/benjbritton/replicalm).
