# =============================================================================
# 05_clean_occurrences.R
# Cleaning of the occurrences and input data of the models.
#
# Takes data/syntrichia_occ_with_climate_vif.csv (03_vif_selection.R) and applies:
#   0b. Exclusion of Syntrichiadelphus flagellaris (outside Syntrichia, Brinda et al. 2021)
#   1. Confidence level high + medium
#   2. Larrain records without coordinates (they came out as 0,0): replaced by the
#      georeference in data/larrain_georeferences.csv if the uncertainty radius is
#      <= 5 km. The rest are excluded.
#   3. Remaining invalid coordinates (0,0 / lat == lon / NA) and records without climate
#   4. basisOfRecord: FOSSIL_SPECIMEN and LIVING_SPECIMEN are excluded
#   5. CoordinateCleaner (Zizka et al. 2019) on the GBIF records:
#      capitals, centroids, lat == lon, GBIF headquarters, institutions, zeros
#   6. Exact duplicates (same taxon and same coordinates)
#   7. Species-level models: only records inside the study area (with values in
#      stack_final_chile.tif = Chile + 200 km)
#   8. Thinning: 1 record per 30-arcsec cell (CHELSA resolution), by species in the
#      species-level models and over the whole genus for the genus-level model
#   S. anderssonii is also saved with its whole range (occ_anderssonii_full_range.csv)
#   to calibrate without truncating the niche (Chevalier et al. 2021).
#
# Inputs:  data/syntrichia_occ_with_climate_vif.csv, data/larrain_georeferences.csv,
#          data/chile_layers/stack_final_chile.tif (04_chile_layers.R)
# Outputs in data/clean/:
#   occ_genus.csv, occ_<species>.csv, occ_anderssonii_full_range.csv (with recordID)
#   cleaning_table.csv     (n per step)
#   excluded_records.csv   (recordID, step and reason of every exclusion, used by
#                           06_figshare_dataset.R)
# Requires: data.table, terra, CoordinateCleaner
# =============================================================================

library(data.table)
library(terra)

# ----------------------------- SETTINGS -------------------------------------
path_csv      <- "data/syntrichia_occ_with_climate_vif.csv"
path_georef   <- "data/larrain_georeferences.csv"
path_stack    <- "data/chile_layers/stack_final_chile.tif"   # 7 variables, Chile + 200 km
out_dir       <- "data/clean"
max_radius_km <- 5
# The 5 species with most records in Chile after cleaning
sp_names      <- c(epilosa       = "Syntrichia epilosa",
                   anderssonii   = "Syntrichia anderssonii",
                   pseudorobusta = "Syntrichia pseudorobusta",
                   breviseta     = "Syntrichia breviseta",
                   squarripila   = "Syntrichia squarripila")
# Species whose range lies mostly outside Chile. Besides the normal file (Chile + 200 km),
# occ_<species>_full_range.csv is saved with the whole range (no area filter), to
# calibrate without truncating the niche (Chevalier et al. 2021).
# S. anderssonii: 217 of 360 records fall outside the area (New Zealand, Australia,
# subantarctic islands).
full_range_species <- c(anderssonii = "Syntrichia anderssonii")
# Taxa outside Syntrichia (Brinda et al. 2021): excluded from the whole analysis
outside_genus <- c("Syntrichia flagellaris")   # = Syntrichiadelphus flagellaris
vif_vars      <- c("bio02", "bio04", "bio08", "bio09", "bio15", "bio18", "bio19")
# ------------------------------------------------------------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
occ <- fread(path_csv, encoding = "UTF-8")
occ[, original_id := trimws(as.character(original_id))]

tab      <- list()
excluded <- list()
register <- function(step, d) {
  row_t <- data.table(step = step, genus = nrow(d))
  for (nm in names(sp_names)) row_t[[nm]] <- d[model_taxon == sp_names[[nm]], .N]
  tab[[length(tab) + 1]] <<- row_t
  cat(sprintf("%-55s %7d records\n", step, nrow(d)))
}
annotate_excluded <- function(x, step, reason)
  excluded[[length(excluded) + 1]] <<- x[, .(recordID, source, original_id, model_taxon,
    decimalLatitude, decimalLongitude, countryCode, step = step, reason = reason)]
exclude <- function(d, drop_rows, step, reason) {
  if (any(drop_rows)) annotate_excluded(d[drop_rows], step, reason)
  d[!drop_rows]
}

register("0. Input file", occ)

# --- 0b. Syntrichiadelphus flagellaris (Brinda et al. 2021) outside the genus ---
# standard_name is used so that doubtful determinations (cf., ?) are also removed
outside <- occ$standard_name %in% outside_genus
cat("Syntrichiadelphus flagellaris records excluded:", sum(outside), "\n")
occ <- exclude(occ, outside, "0b", "Syntrichiadelphus flagellaris (outside Syntrichia, Brinda et al. 2021)")
register("0b. Without Syntrichiadelphus flagellaris", occ)

# --- 1. Confidence high + medium ---------------------------------------------
d <- occ[confidence_level %in% c("high", "medium")]
register("1. Confidence high + medium", d)

# --- 2. Georeferences for Larrain records ------------------------------------------
no_coords <- d$source == "Larrain_herbarium" &
             d$decimalLatitude == 0 & d$decimalLongitude == 0
cat("Larrain records without coordinates:", sum(no_coords), "\n")

n_georef_ok <- 0L
if (file.exists(path_georef)) {
  geo <- fread(path_georef, encoding = "UTF-8")
  geo[, original_id := trimws(as.character(original_id))]
  geo <- geo[!is.na(lat_dec) & !is.na(lon_dec) &
             !is.na(uncertainty_radius_km) & uncertainty_radius_km <= max_radius_km]
  if (nrow(geo) > 0) {
    stk <- rast(path_stack)[[vif_vars]]
    for (i in seq_len(nrow(geo))) {
      idx <- which(no_coords & d$original_id == geo$original_id[i] &
                   d$original_name == geo$original_name[i])
      if (length(idx) == 0) next
      d[idx, `:=`(decimalLatitude  = geo$lat_dec[i],
                  decimalLongitude = geo$lon_dec[i],
                  coordinateUncertaintyInMeters = geo$uncertainty_radius_km[i] * 1000)]
      # Climate re-extracted at the new coordinate (same scale as the CSV, because
      # stack_final_chile.tif went through the same transformation)
      env <- terra::extract(stk, cbind(geo$lon_dec[i], geo$lat_dec[i]))
      for (v in vif_vars) set(d, idx, v, env[[v]])
      no_coords[idx] <- FALSE
      n_georef_ok <- n_georef_ok + length(idx)
    }
  }
} else {
  cat("[WARNING]", path_georef, "does not exist - records without coordinates are excluded.\n")
}
cat("Georeferenced and included (radius <=", max_radius_km, "km):", n_georef_ok, "\n")
d <- exclude(d, no_coords, "2", "Larrain without coordinates, no georeference <= 5 km")
register("2. Larrain without coordinates (georef. or excluded)", d)

# --- 3. Invalid coordinates -----------------------------------------------------------
bad <- is.na(d$decimalLatitude) | is.na(d$decimalLongitude) |
       (d$decimalLatitude == 0 & d$decimalLongitude == 0) |
       (d$decimalLatitude == d$decimalLongitude) |
       abs(d$decimalLatitude) > 90 | abs(d$decimalLongitude) > 180
d <- exclude(d, bad, "3", "Invalid coordinate (0,0 / lat==lon / NA)")
no_climate <- !complete.cases(d[, ..vif_vars])
d <- exclude(d, no_climate, "3", "No climate values")
register("3. Invalid coordinates and no climate", d)

# --- 4. basisOfRecord ----------------------------------------------------------------
bor <- d$basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")
d <- exclude(d, bor, "4", "basisOfRecord fossil or living")
register("4. Exclude fossil and living specimens", d)

# --- 5. CoordinateCleaner (GBIF only, Larrain is expert-curated) ---------------------
if (!requireNamespace("CoordinateCleaner", quietly = TRUE)) {
  if (!identical(Sys.getenv("ALLOW_NO_CC"), "1"))
    stop("CoordinateCleaner is missing: install.packages('CoordinateCleaner')")
  cat("[WARNING] Running WITHOUT CoordinateCleaner (for tests only).\n")
} else {
  is_gbif <- d$source == "GBIF"
  cc <- CoordinateCleaner::clean_coordinates(
    x = as.data.frame(d[is_gbif]),
    lon = "decimalLongitude", lat = "decimalLatitude",
    species = "model_taxon",
    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"),
    capitals_rad = 10000, centroids_rad = 1000, centroids_detail = "both",
    inst_rad = 100, value = "spatialvalid", verbose = FALSE)
  flag <- rep(FALSE, nrow(d))
  flag[which(is_gbif)] <- !cc$.summary
  d <- exclude(d, flag, "5", "CoordinateCleaner (capital/centroid/institution/GBIF)")
}
register("5. CoordinateCleaner (GBIF)", d)

# --- 6. Exact duplicates ----------------------------------------------------------------
# Larrain records are prioritized (expert curation), then the best coordinate quality.
# The comparison uses model_taxon, so after the S. prostrata = S. serrulata synonymy
# two GBIF records of S. serrulata become duplicates of Larrain specimens and are removed.
quality_priority <- c(high = 1, medium = 2, not_reported = 3, low = 4)
d[, priority := fifelse(source == "Larrain_herbarium", 0, quality_priority[coordinate_quality])]
d[is.na(priority), priority := 5]
setorder(d, priority)
d[, lat4 := round(decimalLatitude, 4)][, lon4 := round(decimalLongitude, 4)]
dup <- duplicated(d, by = c("model_taxon", "lat4", "lon4"))
d <- exclude(d, dup, "6", "Exact duplicate (taxon + coordinate)")
register("6. Exact duplicates", d)

# --- 7. 30 arc-second cell ----------------------------------------------------------------
d[, cell := paste(floor((decimalLatitude + 90) * 120), floor((decimalLongitude + 180) * 120))]

save_set <- function(x, slug) {
  out <- x[, c("decimalLongitude", "decimalLatitude", vif_vars,
               "model_taxon", "source", "original_id", "countryCode", "recordID"), with = FALSE]
  setnames(out, c("decimalLongitude", "decimalLatitude"), c("x", "y"))
  fwrite(out, file.path(out_dir, paste0("occ_", slug, ".csv")))
}

# GENUS: whole genus, 1 record per cell
d_genus <- d[!duplicated(d, by = "cell")]
annotate_excluded(d[duplicated(d, by = "cell")], "8_genus", "Thinning, genus-level model (another record already in the cell)")
save_set(d_genus, "genus")

# Species: study area + 1 record per cell
stk <- rast(path_stack)[[vif_vars]]
row_sp <- list(step = c("7. Study area (Chile + 200 km)", "8. Thinning 1 per 30\" cell"))
res_sp <- list()
for (nm in names(sp_names)) {
  s <- d[model_taxon == sp_names[[nm]]]
  inside <- complete.cases(terra::extract(stk, as.matrix(s[, .(decimalLongitude, decimalLatitude)])))
  annotate_excluded(s[!inside], "7_species", "Outside the study area (disjunct)")
  s <- s[inside]
  n_area <- nrow(s)
  dup_cell <- duplicated(s, by = "cell")
  annotate_excluded(s[dup_cell], "8_species", "Thinning, species-level (cell already represented)")
  s <- s[!dup_cell]
  save_set(s, nm)
  res_sp[[nm]] <- c(n_area, nrow(s))
}
# Species with full range (no area filter), 1 record per cell
for (nm in names(full_range_species)) {
  s <- d[model_taxon == full_range_species[[nm]]]
  s <- s[!duplicated(s, by = "cell")]
  save_set(s, paste0(nm, "_full_range"))
  res_sp[[paste0(nm, "_full_range")]] <- c(NA_integer_, nrow(s))
  cat("Full range", nm, ":", nrow(s), "records (1 per cell)\n")
}
for (k in 1:2) {
  row_t <- data.table(step = row_sp$step[k], genus = if (k == 1) NA_integer_ else nrow(d_genus))
  for (nm in names(res_sp)) row_t[[nm]] <- res_sp[[nm]][k]
  tab[[length(tab) + 1]] <- row_t
}

tab <- rbindlist(tab, fill = TRUE)
fwrite(tab, file.path(out_dir, "cleaning_table.csv"))
fwrite(rbindlist(excluded, fill = TRUE), file.path(out_dir, "excluded_records.csv"))

cat("\n================ Cleaning table ================\n")
print(tab)
cat("\nFiles in", out_dir, "\n")
