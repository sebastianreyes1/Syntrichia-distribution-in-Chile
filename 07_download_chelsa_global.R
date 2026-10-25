# =============================================================================
# 07_download_chelsa_global.R
#
# Only needed for the genus-level model (global background) and for the
# S. anderssonii full-range model.
# Downloads to data/chelsa_global/ the 7 CHELSA 2.1 (1981-2010) layers selected
# by VIF. Run ONCE on a cluster node with internet access (usually the login node).
# About 1-3 GB in total.
#
# If you already have the 19 global GeoTIFFs from 02_extract_climate.R
# (data/chelsa_global/CHELSA_bioXX.tif), you do not need to run this: script 08
# looks for exactly those file names.
#
# IMPORTANT (scale): the files are saved as they come from CHELSA. Script 08
# applies to them the SAME transformation used when the climate was extracted at
# the occurrences (temperatures / 10), so that presences, background and
# projection layers are on the same scale.
# =============================================================================

options(timeout = 1800)
dir.create("data/chelsa_global", showWarnings = FALSE, recursive = TRUE)

vars <- c(2, 4, 8, 9, 15, 18, 19)
base_url <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio%d_1981-2010_V.2.1.tif"

for (v in vars) {
  destination <- sprintf("data/chelsa_global/CHELSA_bio%02d.tif", v)
  url      <- sprintf(base_url, v)
  expected <- tryCatch({
    h  <- curlGetHeaders(url)
    cl <- grep("^content-length", h, ignore.case = TRUE, value = TRUE)
    as.numeric(trimws(sub(".*:", "", tail(cl, 1))))
  }, error = function(e) NA)
  if (file.exists(destination) && !is.na(expected) && file.size(destination) == expected) {
    cat("Already exists:", destination, "\n"); next
  }
  for (attempt in 1:3) {
    cat("Downloading bio", v, "(attempt", attempt, ")...\n")
    ok <- tryCatch({ download.file(url, destination, mode = "wb", method = "libcurl", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (ok && (is.na(expected) || file.size(destination) == expected)) break
  }
  if (!file.exists(destination)) stop("Could not download ", url)
}
cat("Done. Layers in data/chelsa_global/\n")
