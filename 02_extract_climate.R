# =============================================================================
# 02_extract_climate.R
# Downloads the 19 CHELSA 2.1 bioclimatic variables (1981-2010) and extracts their
# values at every record with coordinates. Temperatures are divided by 10 (see the
# note on units in the README).
#
# Input:  data/syntrichia_merged.csv (01_merge_sources.R)
# Output: data/syntrichia_occ_with_climate.csv (175,983 records) and data/chelsa_global/
# Requires: terra >= 1.7, data.table. The 19 global layers (2-6 GB) are downloaded
# once and the extraction is done locally.
# =============================================================================
library(terra)
library(data.table)

# --- 1. Load occurrences ---
occ <- fread("data/syntrichia_merged.csv")
occ <- occ[!is.na(decimalLatitude) & !is.na(decimalLongitude)]

pts <- vect(occ, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

# --- 2. Download the 19 global CHELSA GeoTIFFs to a local directory ---
# Some variables (mostly precipitation, e.g. bio12) are much larger than the
# temperature ones (they can exceed 600 MB), so with R's default timeout (60 s) an
# ordinary connection cannot finish. We raise the timeout, force the libcurl method
# (more reliable than the default one in RStudio) and compare the downloaded size
# with the size reported by the server, retrying if they do not match.
options(timeout = 1800)  # 30 minutes per file, enough even for the largest one

dir.create("data/chelsa_global", showWarnings = FALSE, recursive = TRUE)

base_url <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio%d_1981-2010_V.2.1.tif"
local_paths <- sprintf("data/chelsa_global/CHELSA_bio%02d.tif", 1:19)

tif_is_valid <- function(path) {
  # The size can match the expected one and the TIFF can still be corrupt internally
  # (e.g. if the connection dropped while a block was being written). We confirm it
  # by actually reading a few pixels.
  tryCatch({
    r <- rast(path)
    invisible(spatSample(r, size = 5, method = "random", values = TRUE, na.rm = FALSE))
    TRUE
  }, error = function(e) FALSE)
}

download_with_retries <- function(url, destfile, max_attempts = 5) {
  for (attempt in 1:max_attempts) {
    expected_size <- tryCatch({
      h <- curlGetHeaders(url)
      cl <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
      if (length(cl)) as.numeric(gsub("\\D", "", cl[1])) else NA
    }, error = function(e) NA)

    ok <- tryCatch({
      download.file(url, destfile, mode = "wb", quiet = FALSE, method = "libcurl")
      TRUE
    }, error = function(e) {
      cat("  attempt", attempt, "failed:", conditionMessage(e), "\n")
      FALSE
    })

    if (ok) {
      actual_size <- file.size(destfile)
      size_ok <- is.na(expected_size) || actual_size == expected_size
      if (size_ok && tif_is_valid(destfile)) {
        return(TRUE)
      }
      if (!size_ok) {
        cat("  incomplete size (", actual_size, "of", expected_size, "), retrying...\n")
      } else {
        cat("  correct size but corrupt TIFF, retrying...\n")
      }
    }
    if (file.exists(destfile)) file.remove(destfile)
  }
  stop("Could not download ", url, " after ", max_attempts, " attempts")
}

for (i in 1:19) {
  if (file.exists(local_paths[i])) {
    cat("bio", i, "was already downloaded, skipping\n")
    next
  }
  cat("Downloading bio", i, "...\n")
  download_with_retries(sprintf(base_url, i), local_paths[i])
}

# --- 3. Load the local stack and extract (this is local now, it should take
#         seconds to a couple of minutes, not half an hour) ---
bio_stack <- rast(local_paths)
names(bio_stack) <- sprintf("bio%02d", 1:19)

cat("Extracting bioclimatic variables at", nrow(occ), "points...\n")
env_occ <- extract(bio_stack, pts, ID = FALSE)
setDT(env_occ)  # extract() returns a base data.frame, converted to use data.table syntax

# CHELSA distributes temperature in tenths of a degree Celsius: bio1, 2, 4-11 must be
# divided by 10 to get degrees Celsius (check the metadata if you use a different
# version). bio3 and bio12-19 (isothermality, precipitation) already come in their
# final units.
temp_cols <- c("bio01","bio02","bio04","bio05","bio06","bio07","bio08","bio09","bio10","bio11")
env_occ[, (temp_cols) := lapply(.SD, function(x) x / 10), .SDcols = temp_cols]

occ_env <- cbind(occ, env_occ)

fwrite(occ_env, "data/syntrichia_occ_with_climate.csv")
cat("Done: data/syntrichia_occ_with_climate.csv (", nrow(occ_env), "rows,", ncol(occ_env), "columns )\n")
