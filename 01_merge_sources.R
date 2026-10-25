# =============================================================================
# 01_merge_sources.R
# Merges the GBIF download (already reduced by 00_reduce_gbif.R) with the database
# of herbarium specimens of J. Larrain, and assigns to every record:
#   - a standardized name (standard_name) and the taxon used in the analyses
#     (model_taxon; "Syntrichia sp." when the identification is uncertain)
#   - a confidence level (high / medium / low)
#   - a coordinate quality class based on the reported uncertainty
#   - a flag for possible duplicates of Larrain specimens
# GBIF records that are digitized copies of specimens already in the Larrain
# database (duplicates held at NY and UCh) are removed here (33 records).
#
# Taxonomic decision: S. prostrata is treated as a synonym of S. serrulata
# (Gallego et al. 2009, Muller 2009, GBIF backbone). The original name is kept in
# standard_name, and the merged species is used in model_taxon. As a consequence,
# two GBIF records of S. serrulata become exact duplicates of Larrain specimens
# and are removed at the duplicate step of 05_clean_occurrences.R.
#
# Inputs:  data/occurrence_reduced.csv              (00_reduce_gbif.R)
#          data/Syntrichia_Juan_Larrain_2020.xls    (database of J. Larrain)
# Outputs: data/syntrichia_merged.csv               (175,986 records, with recordID)
#          data/gbif_larrain_duplicates.csv         (33 GBIF records removed, for Data S1)
# Requires: data.table, readxl
#
# NOTE: the column names of the Larrain spreadsheet are in Spanish (Especie, Colector,
# Localidad, Habitat, Altitud, Observaciones, Fecha, Mes, Dia, Det., Herbario ...)
# and the months are Spanish abbreviations (Ene, Abr, Ago, Dic). They are read as
# they are in the source file.
# =============================================================================
library(data.table)

path_gbif    <- "data/occurrence_reduced.csv"
path_larrain <- "data/Syntrichia_Juan_Larrain_2020.xls"
path_output  <- "data/syntrichia_merged.csv"

# ---- 1. GBIF -----------------------------------------------------------------
g <- fread(path_gbif, encoding = "UTF-8", colClasses = list(character = "gbifID"), strip.white = FALSE)
blank <- function(x) { x <- as.character(x); x[!is.na(x) & trimws(x) == ""] <- NA; x }
g[, species := blank(species)]
g[, identificationQualifier := blank(identificationQualifier)]

# Qualifiers that indicate an uncertain identification
uncertain <- function(q) {
  q <- ifelse(is.na(q), "", q)
  grepl("^\\?", q) |
    grepl("\\bc\\.?f\\.?r?\\b", q) |     # cf., cf, c.f., cfr. (not "C.fr." = with sporophytes)
    grepl("(?i)aff|dubious|unconfirmed|uncertain|unknown|^sp\\.$|s\\. lat|several|approaching|subulata", q, perl = TRUE)
}
is_inat <- grepl("iNaturalist", g$datasetName)
unc <- as.numeric(g$coordinateUncertaintyInMeters)
specimen <- c("PRESERVED_SPECIMEN", "MATERIAL_SAMPLE", "LIVING_SPECIMEN", "FOSSIL_SPECIMEN")
geo <- as.logical(g$hasGeospatialIssues)

gbif <- g[, .(
  source = "GBIF", original_id = gbifID,
  original_name = scientificName, standard_name = species,
  id_qualifier = identificationQualifier,
  decimalLatitude, decimalLongitude, elevation = as.character(elevation), countryCode, stateProvince,
  locality, habitat = NA_character_, year, month, day, eventDate = as.character(eventDate),
  recordedBy, identifiedBy, institutionCode_x = institutionCode, basisOfRecord,
  # medium: specimen without geospatial issues and with uncertainty <= 5 km or not reported
  confidence_level = fifelse(basisOfRecord %in% specimen & !geo & (is.na(unc) | unc <= 5000), "medium", "low"),
  occurrenceRemarks = NA_character_,
  # possible duplicate: collected by J. Larrain in Chile (not removed here, see step 3)
  possible_larrain_duplicate = grepl("larra.{1,2}n", recordedBy, ignore.case = TRUE) & countryCode %in% "CL",
  confirmed_larrain_duplicate = FALSE,
  taxonRank, coordinateUncertaintyInMeters, coordinatePrecision, hasGeospatialIssues = geo, issue,
  original_gbif_qualifier = identificationQualifier,
  # taxon used in the analyses: genus level if the identification is uncertain, if GBIF
  # does not assign it to a species, or if it is an iNaturalist observation of the
  # S. ruralis complex
  model_taxon = fifelse(is.na(species) | uncertain(identificationQualifier) |
                          (is_inat & species %in% "Syntrichia ruralis"), "Syntrichia sp.", species),
  coordinate_quality = fcase(is.na(unc), "not_reported", unc <= 1000, "high", unc <= 5000, "medium", default = "low"),
  gbif_geospatial_issue = geo, institutionCode_y = institutionCode, datasetName,
  is_inaturalist = is_inat)]

# ---- 2. Larrain database -------------------------------------------------------
L <- if (grepl("\\.csv$", path_larrain)) fread(path_larrain, encoding = "UTF-8", colClasses = "character") else
  as.data.table(readxl::read_excel(path_larrain, col_types = "text", trim_ws = FALSE))
num <- function(x) suppressWarnings(as.numeric(x))
# columns with accents or symbols are found by pattern, so the code does not depend on the encoding
col <- function(pattern) L[[grep(pattern, names(L))[1]]]
txt <- function(x) blank(x)
name <- txt(L$Especie)
# qualifiers written inside the name
qual <- mapply(function(n) {
  if (is.na(n)) return(NA_character_)
  k <- c(if (grepl("\\?", n)) "?", if (grepl("\\bcf\\.", n)) "cf", if (grepl("\\baff\\.", n)) "aff", if (grepl("\\bsp\\.", n)) "sp")
  if (length(k)) paste(k, collapse = ";") else NA_character_
}, name, USE.NAMES = FALSE)
standard <- trimws(gsub("\\s+", " ", gsub("\\?|\\bcf\\. |\\baff\\. ", "", name)))
standard <- sub("andersonii$", "anderssonii", standard)
months <- c(Ene = 1, Feb = 2, Mar = 3, Abr = 4, May = 5, Jun = 6, Jul = 7, Ago = 8, Sep = 9, Oct = 10, Nov = 11, Dic = 12)  # Spanish month abbreviations used in the source file
iso2 <- function(p) { p <- toupper(p); fcase(grepl("^CH", p), "CL", grepl("^ARG", p), "AR", grepl("^ESP", p), "ES", grepl("^PER", p), "PE", default = NA_character_) }
collector <- txt(L$Colector); with_collector <- txt(L$con)   # 'con' = additional collector

lar <- data.table(
  source = "Larrain_herbarium", original_id = txt(col("^N.{1,2}$")),
  original_name = name, standard_name = standard, id_qualifier = qual,
  decimalLatitude = num(L[["Lat (dec)"]]), decimalLongitude = num(L[["Long (dec)"]]),
  elevation = txt(L$Altitud), countryCode = iso2(txt(col("^Pa.{1,2}s$"))),
  stateProvince = txt(L[["Localidad (1st order)"]]), locality = txt(L$Localidad), habitat = txt(L$Habitat),
  year = num(col("^A.{1,2}o$")), month = unname(months[trimws(txt(L$Mes))]), day = num(L$Dia), eventDate = txt(L$Fecha),
  recordedBy = fifelse(is.na(with_collector), collector, paste(collector, with_collector, sep = ";")),
  identifiedBy = txt(L[["Det."]]), institutionCode_x = txt(L$Herbario), basisOfRecord = "PRESERVED_SPECIMEN",
  confidence_level = "high", occurrenceRemarks = txt(L$Observaciones),
  possible_larrain_duplicate = FALSE, confirmed_larrain_duplicate = FALSE, taxonRank = NA_character_,
  coordinateUncertaintyInMeters = NA_real_, coordinatePrecision = NA_real_, hasGeospatialIssues = NA,
  issue = NA_character_, original_gbif_qualifier = NA_character_,
  model_taxon = fifelse(is.na(qual), standard, "Syntrichia sp."),
  coordinate_quality = "not_reported", gbif_geospatial_issue = FALSE,
  institutionCode_y = NA_character_, datasetName = NA_character_, is_inaturalist = FALSE)

# ---- 2b. Synonymy: S. prostrata = S. serrulata ------------------------------------
# Applied to the taxon used in the analyses only. standard_name keeps the original name.
synonyms <- c("Syntrichia prostrata" = "Syntrichia serrulata")
lar[model_taxon %in% names(synonyms),  model_taxon := unname(synonyms[model_taxon])]
gbif[model_taxon %in% names(synonyms), model_taxon := unname(synonyms[model_taxon])]

# ---- 3. Confirmed duplicates -----------------------------------------------------
# GBIF record collected by Larrain that matches a specimen of his database:
# same standardized taxon, same year and within 1 km. The version in the Larrain
# database is kept (expert curation).
cand <- which(grepl("larra.{1,2}n", gbif$recordedBy, ignore.case = TRUE))
dist_km <- function(la1, lo1, la2, lo2) sqrt(((la2 - la1) * 111)^2 + ((lo2 - lo1) * 111 * cos(la1 * pi / 180))^2)
dup <- vapply(cand, function(i) {
  r <- gbif[i]; s <- lar[standard_name %in% r$standard_name & year %in% r$year & !is.na(decimalLatitude)]
  nrow(s) > 0 && min(dist_km(r$decimalLatitude, r$decimalLongitude, s$decimalLatitude, s$decimalLongitude)) <= 1
}, logical(1))
dup_idx <- cand[dup]
cat("Confirmed duplicates of Larrain specimens removed:", length(dup_idx), "\n")
fwrite(g[cand[dup]], "data/gbif_larrain_duplicates.csv")   # original GBIF rows
if (length(dup_idx)) gbif <- gbif[-dup_idx]

# ---- 4. Merge and save ---------------------------------------------------------------
merged <- rbind(lar, gbif, fill = TRUE)
# recordID: stable identifier of every record, carried through the whole workflow.
# It is the recordID of Data S1 (Figshare)
merged <- cbind(recordID = seq_len(nrow(merged)), merged)
cat("Records in the merged set:", nrow(merged), "(Larrain", nrow(lar), "| GBIF", nrow(gbif), ")\n")
print(merged[, .N, by = confidence_level])
fwrite(merged, path_output)
