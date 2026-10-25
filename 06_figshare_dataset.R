# =============================================================================
# 06_figshare_dataset.R
# Builds the dataset published in Figshare (Data S1): all records of the merged set,
# plus the 33 GBIF records removed because they are copies of specimens of the
# Larrain database, with their cleaning flags and their use in each analysis.
# It also builds Table S1 (records per cleaning step) and the number of records
# per species used in the maps of Chile.
#
# Inputs: data/syntrichia_merged.csv          (01_merge_sources.R)
#         data/gbif_larrain_duplicates.csv    (01_merge_sources.R)
#         data/occurrence_reduced.csv         (00_reduce_gbif.R, extra GBIF columns)
#         data/clean/excluded_records.csv and occ_*.csv (05_clean_occurrences.R)
#         data/natural_earth/ne_10m_admin_0_countries.shp (Natural Earth 1:10m)
# Outputs in data/figshare/:
#   Syntrichia_occurrences_with_cleaning_flags.csv  (176,019 records, 50 columns)
#   cleaning_steps_record_counts.csv                 (Table S1)
#   Chile_records_per_species.csv
# Requires: data.table, readr, sf
# =============================================================================
library(data.table)
library(sf)

dir_clean <- "data/clean"
dir_out   <- "data/figshare"
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
model_species <- c(Sepilosa = "Syntrichia epilosa", Sanderssonii = "Syntrichia anderssonii",
                  Spseudorobusta = "Syntrichia pseudorobusta", Sbreviseta = "Syntrichia breviseta",
                  Ssquarripila = "Syntrichia squarripila")
slug <- c(Sepilosa = "epilosa", Sanderssonii = "anderssonii", Spseudorobusta = "pseudorobusta",
          Sbreviseta = "breviseta", Ssquarripila = "squarripila")

# readr keeps trailing spaces and inner quotes of the text fields as they come in
# the sources
read_all <- function(f) as.data.table(readr::read_csv(f, col_types = readr::cols(.default = "c"),
                                                      trim_ws = FALSE, na = "", progress = FALSE))
fus  <- read_all("data/syntrichia_merged.csv")
gb   <- read_all("data/occurrence_reduced.csv")
dupl <- read_all("data/gbif_larrain_duplicates.csv")
exc  <- read_all(file.path(dir_clean, "excluded_records.csv"))
exc[, recordID := as.integer(recordID)]
yes_no <- function(x) toupper(x) %in% c("TRUE", "T")

# ---- 1. Descriptive columns -------------------------------------------------------
is_gbif <- fus$source == "GBIF"
extra <- gb[match(fus$original_id, gbifID), .(collectionCode, catalogNumber, occurrenceStatus)]
extra[!is_gbif, names(extra) := NA]
rename_taxon <- function(x) sub("^Syntrichia flagellaris$", "Syntrichiadelphus flagellaris", x)
binomial <- function(x) {    # analysis species (varieties are assigned to their species)
  w <- strsplit(x, " ", fixed = TRUE)
  vapply(w, function(p) {
    if (length(p) < 2 || is.na(p[2]) || p[2] == "sp.") return(NA_character_)
    if (length(p) > 2 && !p[3] %in% c("var.", "subsp.", "f.")) return(NA_character_)
    paste(p[1:2], collapse = " ")
  }, "")
}
quality <- c(high = "high (<= 1 km)", medium = "medium (1-5 km)", low = "low (> 5 km)", not_reported = "not reported")
confidence <- c(high = "high", medium = "medium", low = "low")

ds <- fus[, .(
  recordID = as.integer(recordID),
  source = fifelse(is_gbif, "GBIF", "Larrain database"),
  sourceRecordID = original_id, gbifID = fifelse(is_gbif, original_id, NA_character_),
  verbatimScientificName = original_name, standardizedName = rename_taxon(standard_name),
  identificationQualifier = id_qualifier, taxonRank,
  analysisTaxon = rename_taxon(model_taxon), analysisSpecies = NA_character_,
  decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, coordinatePrecision,
  coordinateQuality = unname(quality[coordinate_quality]),
  countryCode, stateProvince, locality, habitat, elevation, eventDate, year, month, day,
  recordedBy, identifiedBy, institutionCode = institutionCode_x,
  collectionCode = extra$collectionCode, catalogNumber = extra$catalogNumber, datasetName,
  basisOfRecord, occurrenceStatus = extra$occurrenceStatus, occurrenceRemarks,
  isINaturalist = yes_no(is_inaturalist), gbifHasGeospatialIssues = yes_no(gbif_geospatial_issue),
  gbifIssues = issue, gbifPossibleDuplicateOfLarrain = yes_no(possible_larrain_duplicate),
  confidenceLevel = unname(confidence[confidence_level]))]
ds[, analysisSpecies := binomial(analysisTaxon)]

# GBIF copies of Larrain specimens (removed in 01_merge_sources.R)
dp <- dupl[, .(source = "GBIF", sourceRecordID = gbifID, gbifID,
  verbatimScientificName = scientificName, standardizedName = species,
  identificationQualifier, taxonRank, decimalLatitude, decimalLongitude,
  coordinateUncertaintyInMeters, coordinatePrecision, countryCode, stateProvince, locality,
  elevation, eventDate, year, month, day, recordedBy, identifiedBy, institutionCode,
  collectionCode, catalogNumber, datasetName, basisOfRecord, occurrenceStatus, gbifIssues = issue)]
dp[, recordID := nrow(ds) + seq_len(.N)]
ds <- rbind(ds, dp, fill = TRUE)

# ---- 2. Exclusion step ---------------------------------------------------------------
reasons <- c(
  "0"  = "GBIF copy of a specimen already in the Larrain database (removed when merging sources)",
  "0b" = "Syntrichiadelphus flagellaris, not a member of Syntrichia (Brinda et al. 2021)",
  "1"  = "Low confidence level (observation-based record, or specimen with GBIF geospatial issue or coordinate uncertainty > 5 km)",
  "2"  = "Larrain record without coordinates, locality not georeferenceable with uncertainty <= 5 km",
  "3"  = "Invalid coordinates or missing climate values",
  "4"  = "Fossil or living specimen",
  "5"  = "Flagged by CoordinateCleaner (capital, centroid, institution or GBIF HQ)",
  "6"  = "Exact duplicate (same taxon and coordinates to 4 decimals)")
ds[, exclusionStep := NA_character_]
ds[recordID > nrow(fus), exclusionStep := "0"]
step_exc <- exc[step %in% names(reasons)]
ds[step_exc, on = .(recordID = recordID), exclusionStep := i.step]
ds[is.na(exclusionStep) & confidenceLevel == "low", exclusionStep := "1"]
# Larrain records without coordinates in the database (they do not reach 05_clean_occurrences.R)
ds[is.na(exclusionStep) & source == "Larrain database" & is.na(decimalLatitude), exclusionStep := "2"]
ds[, exclusionReason := unname(reasons[exclusionStep])]
ds[, passedCleaning := is.na(exclusionStep)]

# ---- 3. Use in the models ----------------------------------------------------------------
ids <- function(f) as.integer(fread(file.path(dir_clean, f), select = "recordID")$recordID)
ds[, usedSDM_genus := recordID %in% ids("occ_genus.csv")]
for (k in names(model_species)) set(ds, j = paste0("usedSDM_", k), value = ds$recordID %in% ids(paste0("occ_", slug[k], ".csv")))
exc_id <- function(p) exc[step == p, recordID]
ds[, sdmGenusNote := fcase(usedSDM_genus, "used",
  recordID %in% exc_id("8_genus"), "not used: 30-arcsec cell already represented (genus-level thinning)",
  default = NA_character_)]
used_sp <- Reduce(`|`, ds[, paste0("usedSDM_", names(model_species)), with = FALSE])
ds[, sdmSpeciesNote := fcase(used_sp, "used",
  recordID %in% exc_id("7_species"), "not used: outside study area (Chile + 200 km)",
  recordID %in% exc_id("8_species"), "not used: 30-arcsec cell already represented (species-level thinning)",
  default = NA_character_)]

# ---- 4. Maps of Chile ---------------------------------------------------------------------
# Cleaned records with countryCode CL and inside continental Chile (5 km coastal
# tolerance, measured in UTM 19S)
sf_use_s2(FALSE)
chile <- st_read("data/natural_earth/ne_10m_admin_0_countries.shp", quiet = TRUE)
chile <- st_geometry(chile[chile$ADMIN == "Chile", ])
chile <- suppressWarnings(suppressMessages(st_crop(chile, st_bbox(c(xmin = -76, ymin = -56.5, xmax = -66, ymax = -17), crs = st_crs(4326)))))
cand <- ds[passedCleaning & countryCode %in% "CL"]
p <- st_as_sf(cand, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
inside <- lengths(st_is_within_distance(st_transform(p, 32719), st_transform(chile, 32719), 5000)) > 0
ds[, usedChileMaps := recordID %in% cand$recordID[inside]]

column_order <- c("recordID", "source", "sourceRecordID", "gbifID", "verbatimScientificName", "standardizedName",
  "identificationQualifier", "taxonRank", "analysisTaxon", "analysisSpecies", "decimalLatitude",
  "decimalLongitude", "coordinateUncertaintyInMeters", "coordinatePrecision", "coordinateQuality",
  "countryCode", "stateProvince", "locality", "habitat", "elevation", "eventDate", "year", "month", "day",
  "recordedBy", "identifiedBy", "institutionCode", "collectionCode", "catalogNumber", "datasetName",
  "basisOfRecord", "occurrenceStatus", "occurrenceRemarks", "isINaturalist", "gbifHasGeospatialIssues",
  "gbifIssues", "gbifPossibleDuplicateOfLarrain", "confidenceLevel", "passedCleaning", "exclusionStep",
  "exclusionReason", "usedSDM_genus", paste0("usedSDM_", names(model_species)), "usedChileMaps",
  "sdmGenusNote", "sdmSpeciesNote")
setcolorder(ds, column_order)
fwrite(ds, file.path(dir_out, "Syntrichia_occurrences_with_cleaning_flags.csv"))

# ---- 5. Table S1: records per step ------------------------------------------------------------
row_fun <- function(label, x, all = TRUE, chile = TRUE) {
  r <- data.table(step = label, all_records = if (all) nrow(x) else NA_integer_)
  for (k in names(model_species)) set(r, j = sub("^S", "S_", k), value = x[analysisTaxon == model_species[[k]], .N])
  r[, chile_maps := if (chile) x[countryCode %in% "CL", .N] else NA_integer_]
  r
}
steps <- c("0b" = "0b. Syntrichiadelphus flagellaris removed (Brinda et al. 2021)",
           "1" = "1. Confidence level high or medium",
           "2" = "2. Larrain records without coordinates removed (not georeferenceable <= 5 km)",
           "3" = "3. Invalid coordinates or missing climate values removed",
           "4" = "4. Fossil and living specimens removed",
           "5" = "5. CoordinateCleaner flags removed (GBIF records only)",
           "6" = "6. Exact duplicates removed (taxon + coordinates)")
x <- ds[exclusionStep %in% NA | exclusionStep != "0"]
tab <- list(row_fun("0. Merged dataset (GBIF records with coordinates + Larrain database)", x))
for (s in names(steps)) { x <- x[exclusionStep %in% NA | exclusionStep != s]; tab[[length(tab) + 1]] <- row_fun(steps[[s]], x) }
f7 <- row_fun("7. Species models: within study area (Chile + 200 km)", x[!recordID %in% exc_id("7_species")], all = FALSE, chile = FALSE)
f8 <- data.table(step = "8. Spatial thinning, 1 record per 30-arcsec cell (final SDM input)",
                 all_records = ds[, sum(usedSDM_genus)])
for (k in names(model_species)) set(f8, j = sub("^S", "S_", k), value = sum(ds[[paste0("usedSDM_", k)]]))
f8[, chile_maps := NA_integer_]
f9 <- data.table(step = "9. Chile maps: countryCode CL and within continental Chile (5 km coastal tolerance)",
                 chile_maps = ds[, sum(usedChileMaps)])
tab <- rbindlist(c(tab, list(f7, f8, f9)), fill = TRUE)
fwrite(tab, file.path(dir_out, "cleaning_steps_record_counts.csv"))
print(tab)

# ---- 6. Records per species in the maps of Chile -------------------------------------------------
ch <- ds[usedChileMaps == TRUE]
ch[, species := fifelse(is.na(analysisSpecies), "Syntrichia sp. (genus-level or uncertain identification)", analysisSpecies)]
by_sp <- ch[, .(records = .N, larrain = sum(source == "Larrain database"), gbif = sum(source == "GBIF")), by = species]
setorder(by_sp, -records, species)
by_sp <- rbind(by_sp[!grepl("^Syntrichia sp\\.", species)], by_sp[grepl("^Syntrichia sp\\.", species)])
fwrite(by_sp, file.path(dir_out, "Chile_records_per_species.csv"))
cat("\nDataset:", nrow(ds), "records |", ds[, sum(passedCleaning)], "cleaned |",
    ds[, sum(usedChileMaps)], "in the maps of Chile |", nrow(by_sp) - 1, "species\n")
