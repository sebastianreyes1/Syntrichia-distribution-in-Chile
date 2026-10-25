# ============================================================================
# 11_download_protected_area_layers.R
# Downloads the protected-area layers of Chile used in Fig. 4 and Table S3
# ----------------------------------------------------------------------------
# Source: Spatial Data Infrastructure (IDE) of the Superintendencia del Medio
# Ambiente (public ArcGIS REST server), which replicates layers of MMA / SINIA /
# CONAF / Bienes Nacionales.
#   https://ideserver.sma.gob.cl/arcgis/rest/services/IDE
# Each layer is downloaded in chunks of 5 objects (see download_layer) and saved
# as a GeoPackage in `protected_area_layers/`.
#
# Layers used in the manuscript (fig_records_protected_areas.R and make_figs_final.py):
#   A) State protected areas (formerly SNASPE, MBN 2024)  -> state_pa_mbn
#   B) Nature sanctuaries (MMA)                           -> nature_sanctuaries
#   C) Private protected areas (MMA)                      -> private_pa_a + private_pa_b
#   D) Priority sites (MMA 2020)                          -> priority_sites_mma2020
# priority_sites_conama is downloaded but not used.
# Requires: httr, jsonlite, sf
# ============================================================================
library(httr)
library(jsonlite)
library(sf)
sf_use_s2(FALSE)   # avoids validity errors in large polygons

base    <- "https://ideserver.sma.gob.cl/arcgis/rest/services/IDE"
out_dir <- "protected_area_layers"
dir.create(out_dir, showWarnings = FALSE)

layers <- c(
  # a) State protected areas: national parks, national reserves and natural monuments
  #    (Ministerio de Bienes Nacionales 2024, 107 units)
  state_pa_mbn         = "Areas_proteccion_oficial2/MapServer/41",
  # b) Nature sanctuaries (MMA)
  nature_sanctuaries   = "Areas_proteccion_oficial2/MapServer/40",
  # c) Private protected areas (MMA) and private and community conservation areas (MMA)
  private_pa_a         = "Otras_bajo_proteccion/MapServer/2",
  private_pa_b         = "Otras_bajo_proteccion/MapServer/3",
  # d) Priority sites for biodiversity conservation
  priority_sites_mma2020 = "Otras_bajo_proteccion/MapServer/7",   # current MMA version
  priority_sites_conama  = "Biodiversidad/MapServer/3"            # National Biodiversity Strategy (CONAMA)
  # CAUTION: Biodiversidad/MapServer/1 ("Areas de Preservacion Ecologica") is NOT the state
  # system of protected areas: it holds 40 zones of regulatory plans. If a file
  # protected_area_layers/snaspe.gpkg exists from earlier tests, delete it.
)

# ---- download in chunks of objects ----------------------------------------------
# The server is ArcGIS 10.31: it only returns Esri JSON (f = json), which GDAL reads.
# Some layers (e.g. the state protected areas) have very detailed polygons and the
# server fails or takes very long if they are requested whole. Therefore:
#   - the OBJECTIDs are requested first and downloaded in chunks of 5, with a 2 min limit
#   - the geometry is requested simplified to ~10 m (maxAllowableOffset = 0.0001
#     degrees). For overlaying points with uncertainty of up to 5 km this changes nothing
#   - if an object still fails, it is requested simplified to ~100 m (0.001 degrees)
read_esrijson <- function(txt) {
  tmp <- tempfile(fileext = ".json")
  writeLines(txt, tmp, useBytes = TRUE)
  st_read(paste0("ESRIJSON:", tmp), quiet = TRUE)
}
request_chunk <- function(url, ids, offset = 0.0001) {
  q <- list(objectIds = paste(ids, collapse = ","), outFields = "*", outSR = 4326,
            returnGeometry = "true", f = "json")
  if (!is.null(offset)) q$maxAllowableOffset <- offset
  # timeout: if the server does not answer within 2 min, it is cut and retried with a smaller chunk
  r <- POST(url, body = q, encode = "form", timeout(120))
  txt <- content(r, as = "text", encoding = "UTF-8")
  if (http_error(r) || grepl('"error"', substr(txt, 1, 200))) stop(substr(txt, 1, 200))
  read_esrijson(txt)
}
download_layer <- function(path, chunk = 5) {
  url <- paste0(base, "/", path, "/query")
  ids <- content(GET(url, query = list(where = "1=1", returnIdsOnly = "true", f = "json")),
                 as = "parsed", type = "application/json")$objectIds
  ids <- sort(unlist(ids))
  message(path, ": ", length(ids), " objects")
  if (length(ids) == 0) return(NULL)
  parts <- list(); simplified <- 0
  for (g in split(ids, ceiling(seq_along(ids) / chunk))) {
    x <- tryCatch(request_chunk(url, g), error = function(e) NULL)
    if (is.null(x)) {                       # failed chunk: one object at a time, more simplified
      x <- do.call(rbind, lapply(g, function(id) {
        y <- tryCatch(request_chunk(url, id), error = function(e) NULL)
        if (is.null(y)) { y <- request_chunk(url, id, offset = 0.001); simplified <<- simplified + 1 }
        y
      }))
    }
    parts[[length(parts) + 1]] <- x
    message("  ", sum(sapply(parts, nrow)), " / ", length(ids), " objects")
  }
  x <- do.call(rbind, parts)
  if (simplified > 0) message("  ", simplified, " objects downloaded with geometry simplified to ~100 m")
  if (nrow(x) != length(ids)) warning(path, ": expected ", length(ids), " objects and got ", nrow(x))
  st_crs(x) <- 4326
  st_make_valid(x)
}

for (nm in names(layers)) {
  f <- file.path(out_dir, paste0(nm, ".gpkg"))
  if (file.exists(f)) { message("already exists: ", f); next }
  x <- tryCatch(download_layer(layers[[nm]]), error = function(e) { message("FAILED ", nm, ": ", conditionMessage(e)); NULL })
  if (!is.null(x)) {
    st_write(x, f, delete_dsn = TRUE, quiet = TRUE)
    message("saved: ", f, " (", nrow(x), " objects)")
  }
}

# ---- inspection: what each layer contains ----------------------------------------
# Used to decide which rows correspond to each type in the manuscript. The attribute
# names of the layers are in Spanish (tipo, categoria, clasificacion, nombre), which
# is why the search pattern includes Spanish words.
for (f in list.files(out_dir, pattern = "\\.gpkg$", full.names = TRUE)) {
  x <- st_read(f, quiet = TRUE)
  cat("\n=====", basename(f), "|", nrow(x), "objects | fields:\n")
  print(names(x))
  # first text fields that look like a category (type / category / classification / name)
  cand <- grep("tipo|categ|clasif|nombre|name", names(x), ignore.case = TRUE, value = TRUE)
  for (cc in head(cand, 3)) { cat("\n--", cc, "\n"); print(head(sort(table(x[[cc]], useNA = "ifany"), decreasing = TRUE), 15)) }
}
