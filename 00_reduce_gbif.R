# =============================================================================
# 00_reduce_gbif.R
# Reduces occurrence.txt from the GBIF download (DOI 10.15468/dl.v3dfb4) to the
# columns used in the analyses and drops records without coordinates.
#
# Input:  data/GBIF_Darwin_Core_Syntrichia_Global/occurrence.txt
# Output: data/occurrence_reduced.csv (174,705 records)
# Requires: data.table
# =============================================================================
library(data.table)

# 1. Path to occurrence.txt (GBIF download: https://doi.org/10.15468/dl.v3dfb4)
path_occurrence <- "data/GBIF_Darwin_Core_Syntrichia_Global/occurrence.txt"

columns <- c(
  "gbifID", "scientificName", "acceptedScientificName", "taxonRank", "taxonomicStatus",
  "genus", "species", "specificEpithet", "infraspecificEpithet", "identificationQualifier",
  "typeStatus", "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
  "coordinatePrecision", "countryCode", "stateProvince", "locality", "elevation",
  "elevationAccuracy", "basisOfRecord", "individualCount", "eventDate", "year", "month", "day",
  "institutionCode", "collectionCode", "catalogNumber", "datasetName", "recordedBy",
  "identifiedBy", "dateIdentified", "occurrenceStatus", "hasCoordinate", "hasGeospatialIssues",
  "issue"
)

# The file is tab-separated, with a header in the first line
dt <- fread(
  path_occurrence,
  sep = "\t",
  select = columns,
  encoding = "UTF-8",
  quote = "",
  na.strings = c("", "NA")
)

cat("Total rows:", nrow(dt), "\n")

# Drop records without coordinates (they cannot be used in distribution models)
dt_with_coords <- dt[!is.na(decimalLatitude) & !is.na(decimalLongitude)]
cat("Rows with coordinates:", nrow(dt_with_coords), "\n")

output <- "data/occurrence_reduced.csv"
fwrite(dt_with_coords, output)
cat("Saved to:", output, "\n")
