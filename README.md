# MayaBRIC

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23924657.svg)](https://doi.org/10.5281/zenodo.23924657)

## In plain terms

**What this is.** MayaBRIC gives every ancient Maya structure on a map a score from
0 to 1. The score describes how well the place was connected to salt and to the
political centres that organised its trade. This repository holds the code that
computes the score, the settings it uses, and the locations of the twelve salt
sources it measures against.

**Why it exists.** Salt was a daily necessity for the ancient Maya, for diet and for
preserving food. Most of it came from a few places, mainly on the coast, but many
people lived far inland and depended on trade to get it. A household close to salt
and to several large centres had more ways to get salt than one far from both, so
it was less exposed when one route or one source failed. Airborne LiDAR now shows
tens of thousands of buildings that were never on any map. In our study, a
deep-learning model found 57,181 of them on the Yucatán Peninsula. MayaBRIC applies
the same open, repeatable measure to every one of them. Researchers can then compare
places across the whole region instead of one site at a time. Anyone can also check
the published numbers, because this code reproduces them.

**What it does.** For any point on the map, it asks three questions:

1. **How far is the nearest salt?** It works out the travel time, in hours, on foot
   or by canoe from the point to the cheapest of the twelve salt sources. Travel
   accounts for hills, rivers, the coast and water that can't be crossed. (Pillar A)
2. **How many important centres are nearby?** It counts the ranked sites within
   25 km, about a day's round trip on foot, and gives larger centres more weight.
   (Pillar D)
3. **How strongly is the point tied to a major centre?** It measures the strongest
   pull from one nearby centre, which is stronger when the centre is larger and
   closer. (Pillar G)

Each answer becomes a score between 0 and 1. MayaBRIC is the average of the three,
and the result is sorted into four bands, 0 to 3. The code also groups the
structures into site-level clusters. It flags clusters with no recorded site within
10 km, because these may be settlements that have not been recorded.

**Who it's for.** Archaeologists and GIS researchers who want to use or extend the
scores, and reviewers or readers who want to check how the published catalog was
made.

**How to use it.** There are three common uses:

- *Score your own points.* Start with any table of point locations (latitude and
  longitude). Run steps 1–3 under [Use](#use) to get travel time to salt, the three
  pillar scores, the MayaBRIC score and band, and the clusters.
- *Score a single named site.* `scripts/score_site.py` takes a name and coordinates
  and returns that site's scores. For example, Calakmul scores 0.591.
- *Check the published catalog.* `tests/verify_against_catalog.py` recomputes every
  score and compares it with the catalog.

**What you need, and what isn't here.** Some inputs can't be shared in this
repository. The locations of the 57,181 structures are restricted to protect the
sites from looting. Identified researchers can get them from the corresponding
author. The ranked site list comes from the Electronic Atlas of Ancient Maya Sites
and is available from its authors. The travel-cost map is built from free public
data, but the script that builds it is not included yet. Until it is, its settings
are recorded in `config/friction_config.yaml`. [docs/INPUTS.md](docs/INPUTS.md)
lists every input and where to get it.

## Technical summary

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
