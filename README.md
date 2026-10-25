# Code: *Syntrichia* in Chile. Occurrence data, distribution models and protected areas

This repository contains the code to reproduce the analyses, figures and tables of the manuscript *Are mosses protected where they occur? Distribution and sampling gaps in the genus Syntrichia (Pottiaceae: Bryophyta) in Chile* (Reyes, Fontúrbel and Atala). The data and models are built in **R**; the maps of records, richness and protected areas (Figs. 1, 2, 4 and S1) and Table S3 are drawn in **Python**. Run every script from the project root (the folder that contains this README), in the order below.

Data S1 (every record with its cleaning and use flags) is deposited at [Figshare DOI to be added]. Script `06` rebuilds it from the raw sources.

## Quick map: manuscript item → script

| Manuscript item | Script | File written |
|---|---|---|
| Data S1 | `06_figshare_dataset.R` | `data/figshare/Syntrichia_occurrences_with_cleaning_flags.csv` |
| Table S1 (cleaning steps) | `06_figshare_dataset.R` | `data/figshare/cleaning_steps_record_counts.csv` |
| Fig. 1, Fig. 2, Fig. 4, Fig. S1 | `python/make_figures.py` | `figures/Fig1_…`, `Fig2_…`, `Fig4_…`, `FigS1_…` (PNG, PDF, TIFF) |
| Table S3 (protected areas) | `python/make_figures.py` | `figures/TableS3_species_in_protected_areas.csv` |
| Table 1, Table S2 | `14_model_performance_tables.R` | `figures/Table1_modelled_taxa.csv`, `figures/TableS2_model_performance.csv` |
| Fig. 3 and Fig. S2 | `fig_sdm_suitability.R` | `figures/Fig3_SDM_suitability`, `figures/FigS2_SDM_suitability_records` |
| Fig. S3 (MESS) | `fig_sdm_suitability.R` with `MESS=1` | `figures/FigS3_MESS` |
| Fig. S4 (sampling bias) | `12_sampling_bias.R`, then `fig_sampling_bias.R` | `figures/FigS4_sampling_bias` |
| Fig. S5 and Table S4 | `09_variable_importance.R`, then `fig_variables.R` | `figures/FigS5_variable_importance_response`, `figures/TableS4_variable_importance.csv` |
| Table S5 (ODMAP) | not code (protocol in the Supplementary Material) | |
| Table S6 (target-group background) | `08_sdm_model.R <taxon> tgb`, then `13_tgb_sensitivity.R` | `figures/TableS6_tgb_sensitivity.csv` |

## 1. Occurrence data and predictors

| Script | What it does | Main outputs |
|---|---|---|
| `00_reduce_gbif.R` | Reduces the GBIF download (`occurrence.txt`, https://doi.org/10.15468/dl.v3dfb4) to the columns used and drops records without coordinates | `data/occurrence_reduced.csv` (174,705) |
| `01_merge_sources.R` | Merges GBIF with the J. Larraín herbarium database. Standardizes names, applies the synonymy *S. prostrata* = *S. serrulata*, assigns the taxon used in the models, confidence level and coordinate quality. Removes GBIF copies of Larraín specimens (same taxon, year and ≤ 1 km) and assigns a `recordID` to every record | `data/syntrichia_merged.csv` (175,986), `data/gbif_larrain_duplicates.csv` (33) |
| `02_extract_climate.R` | Downloads the 19 CHELSA 2.1 bioclimatic variables (1981–2010) and extracts them at each record | `data/syntrichia_occ_with_climate.csv`, `data/chelsa_global/` |
| `03_vif_selection.R` | Stepwise VIF (`usdm::vifstep`, threshold 5). Keeps bio02, bio04, bio08, bio09, bio15, bio18, bio19 | `data/syntrichia_occ_with_climate_vif.csv` |
| `04_chile_layers.R` | Crops the 7 variables to Chile + 200 km (projection area and study area of the species models) | `data/chile_layers/stack_final_chile.tif` |
| `05_clean_occurrences.R` | Cleaning steps 0b–6 (Table S1), georeferencing of Larraín records from `data/larrain_georeferences.csv`, study area and thinning to one record per 30-arcsec cell | `data/clean/occ_*.csv`, `excluded_records.csv`, `cleaning_table.csv` |
| `06_figshare_dataset.R` | Builds Data S1, Table S1 and the records per species in Chile | `data/figshare/` |

After cleaning, 1,857 records fall in continental Chile (1,595 identified to species, 39 species). The genus-level model uses 21,242 worldwide presences.

## 2. Distribution models

| Script | What it does | Main outputs |
|---|---|---|
| `07_download_chelsa_global.R` | Downloads the 7 global layers needed by the genus and *S. anderssonii* models (not needed if `02` was run on the same machine) | `data/chelsa_global/` |
| `08_sdm_model.R` | Fits GLM, GAM, MaxNet, RF and SVM with flexsdm (spatial block cross-validation), keeps the algorithms with TSS ≥ 0.6 and Boyce > 0, averages them and projects onto Chile | `data/results/<taxon>_<mode>/` |
| `09_variable_importance.R` | Permutation variable importance and partial response curves | `variable_importance.csv`, `response_curves.csv` |
| `10_mess.R` | MESS between continental Chile and the calibration data of each model | `mess/` |

Usage: `Rscript 08_sdm_model.R <taxon> <mode>`.

* Taxa: `genus` (worldwide records), `epilosa`, `anderssonii_full_range` (full range of *S. anderssonii*, used in the manuscript), `pseudorobusta`, `breviseta`, `squarripila`. The taxon `anderssonii` fits the Chile + 200 km set and is only a comparison.
* Modes: `main` (results of the manuscript), `tgb` (target-group background sensitivity analysis), `test` (quick test).
* On a SLURM cluster: `sbatch job_r.sh 08_sdm_model.R epilosa main`. Fix the module name and resources inside `job_r.sh`.
* Fully reproducing the models takes several hours per taxon. All random processes use `set.seed(42)`.

Example for the whole set of models:

```
for t in genus epilosa anderssonii_full_range pseudorobusta breviseta squarripila; do
  Rscript 08_sdm_model.R $t main
done
for t in epilosa anderssonii_full_range pseudorobusta breviseta squarripila; do
  Rscript 08_sdm_model.R $t tgb
done
for t in epilosa anderssonii_full_range pseudorobusta breviseta squarripila genus; do
  Rscript 09_variable_importance.R $t
done
```

## 3. Protected areas, sampling bias and sensitivity

| Script | What it does | Main outputs |
|---|---|---|
| `11_download_protected_area_layers.R` | Downloads the protected-area layers (state protected areas, nature sanctuaries, private protected areas, priority sites) from the public ArcGIS server of the Chilean Ministry of the Environment | `protected_area_layers/*.gpkg` |
| `12_sampling_bias.R` | Sampling gaps (50 km hexagons) and accessibility bias (distance to roads and towns, Kadmon et al. 2004) | `figures/sampling_bias/` |
| `13_tgb_sensitivity.R` | Compares main and target-group background predictions (Spearman correlation, Schoener's D) | `figures/TableS6_tgb_sensitivity.csv` |

## 4. Figures and tables

Run after the steps above:

```
python python/make_figures.py                # Figs. 1, 2, 4, S1 and Table S3
Rscript 14_model_performance_tables.R        # Table 1 and Table S2
Rscript fig_sdm_suitability.R                # Fig. 3, Fig. S2 and one map per taxon
MESS=1 Rscript fig_sdm_suitability.R         # Fig. S3 (after 10_mess.R)
Rscript fig_sampling_bias.R                  # Fig. S4 (after 12_sampling_bias.R)
Rscript fig_variables.R                      # Fig. S5 and Table S4 (after 09_variable_importance.R)
```

R figures are written as TIFF (600 dpi, LZW) and PDF in `figures/`. Python figures are written as PNG, PDF and TIFF. The `TGB=1` option of `fig_sdm_suitability.R` draws the maps of the target-group-background models.

## Inputs not included

| Input | Where to get it |
|---|---|
| GBIF download `data/GBIF_Darwin_Core_Syntrichia_Global/occurrence.txt` | https://doi.org/10.15468/dl.v3dfb4 |
| J. Larraín database `data/Syntrichia_Juan_Larrain_2020.xls` | Available from the authors on request. The column names and month abbreviations are in Spanish (see the note in `01_merge_sources.R`) |
| CHELSA 2.1 bioclimatic layers (1981–2010) | https://chelsa-climate.org (downloaded by `02` and `07`) |
| Protected-area layers | Downloaded by `11` |
| Natural Earth 1:10m: `ne_10m_admin_0_countries`, `ne_10m_roads`, `ne_10m_populated_places` | https://www.naturalearthdata.com, unzip into `data/natural_earth/` (the Python script downloads the countries layer on its own) |

`data/larrain_georeferences.csv` lists the Larraín records that had no coordinates, with the decision taken for each. None could be georeferenced within the 5 km uncertainty limit, so all were excluded.

## Software

R 4.4.0 with flexsdm 1.4.0, maxnet, usdm, terra, sf, data.table, dplyr, readxl, CoordinateCleaner, ggplot2, patchwork, rnaturalearthdata, httr, jsonlite. Python 3 with the packages in `python/requirements.txt` (fonts: Liberation Sans, or edit `FONT`/`font` to use Arial).

## Notes

**Taxonomy:** *Syntrichia prostrata* is treated as a synonym of *S. serrulata* (Gallego et al. 2009; Müller 2009; GBIF backbone). The original name is kept in the `standard_name` column of the merged file and in Data S1, and the merged species is used for the models and maps.

**MaxNet:** `flexsdm::fit_max()` called without the `background` argument fits each cross-validation fold with presences and pseudo-absences, but fits the final model with presences only, which gives an empty model. `08_sdm_model.R` refits the final MaxNet model with all presences and pseudo-absences, using the same settings as the folds. Cross-validation metrics are not affected.

**Units:** in the modelling data, the CHELSA temperature variables (bio02, bio04, bio08, bio09) are divided by 10. `fig_variables.R` converts them back for the response curves.

**Record identifiers:** `recordID` is assigned in `01` and carried through every step, so each record in Data S1 can be traced to the models and maps that used it.

**Two languages:** Figs. 1, 2, 4, S1 and Table S3 are produced with Python (geopandas, matplotlib) and the remaining figures with R (terra, sf, ggplot2). Both use the same projection, base map and colours.

## Citation

[Manuscript reference and repository DOI (Zenodo) to be added.]
