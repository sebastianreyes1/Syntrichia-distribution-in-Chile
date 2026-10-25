# =============================================================================
# 04_chile_layers.R
# Builds the projection stack: the 7 CHELSA variables selected by VIF, cropped to
# Chile + 200 km. This is the area where the models are projected and the study
# area of the species-level models.
#
# Input:  CHELSA 2.1 layers (read remotely, only the blocks of the area are downloaded)
# Output: data/chile_layers/stack_final_chile.tif (7 layers, 30 arc-seconds)
# Note:   topography was prepared but not used in the models, so it is omitted here.
# Requires: terra, rnaturalearth (and rnaturalearthdata), sf
# =============================================================================
library(terra)
library(rnaturalearth)
library(sf)

dir.create("data/chile_layers", showWarnings = FALSE, recursive = TRUE)

# --- 1. Crop area: Chile + 200 km buffer ---
chile <- ne_countries(country = "Chile", scale = "medium", returnclass = "sf")
chile_buf <- st_buffer(st_transform(chile, 3857), 200000)  # buffer in metres (Web Mercator)
chile_buf <- st_transform(chile_buf, 4326)
chile_v <- vect(chile_buf)

# --- 2. CHELSA bioclim, remote crop ---
base_url <- "/vsicurl/https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio%d_1981-2010_V.2.1.tif"
bio_stack <- rast(sprintf(base_url, 1:19))
names(bio_stack) <- sprintf("bio%02d", 1:19)
bio_chile <- crop(bio_stack, chile_v, mask = TRUE)

# Same unit adjustment as in 02_extract_climate.R (temperatures / 10)
temp_idx <- c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11)
for (i in temp_idx) bio_chile[[i]] <- bio_chile[[i]] / 10

# --- 3. The 7 variables selected by VIF (03_vif_selection.R) ---
vif_variables <- c("bio02", "bio04", "bio08", "bio09", "bio15", "bio18", "bio19")
stack_final <- bio_chile[[vif_variables]]

writeRaster(stack_final, "data/chile_layers/stack_final_chile.tif", overwrite = TRUE,
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
cat("Done: data/chile_layers/stack_final_chile.tif |", paste(names(stack_final), collapse = ", "), "\n")
