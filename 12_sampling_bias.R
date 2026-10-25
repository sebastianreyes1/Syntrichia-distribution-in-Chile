# 12_sampling_bias.R
# Sampling bias of Syntrichia in continental Chile
# (1) Sampling gaps: 50 km hexagons with high suitability for the genus, analogous
#     climate (MESS >= 0) and no records.
# (2) Accessibility bias (Kadmon et al. 2004): distance from the collection
#     localities to roads and populated places vs. random points sampled in
#     proportion to the suitability of the genus (what is expected if collecting
#     followed climate only).
# Inputs:  data/figshare/Syntrichia_occurrences_with_cleaning_flags.csv (06_figshare_dataset.R)
#          data/results/genus_main/ensemble_mean.tif and model_metrics.csv (08_sdm_model.R)
#          mess/mess_genus.tif (10_mess.R)
#          data/natural_earth/ (ne_10m_admin_0_countries, ne_10m_roads, ne_10m_populated_places,
#          from naturalearthdata.com)
# Outputs: figures/sampling_bias/ (CSV tables and bias.rds, used by fig_sampling_bias.R)
suppressPackageStartupMessages({ library(sf); library(terra); library(dplyr) })
sf_use_s2(FALSE); set.seed(42)
albers <- "+proj=aea +lat_1=-22 +lat_2=-50 +lat_0=-36 +lon_0=-71 +datum=WGS84 +units=m"
R <- "."   # project folder
NE <- "data/natural_earth"
OUT <- "figures/sampling_bias"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
thr_rf <- as.numeric(read.csv(file.path(R, "data/results/genus_main/model_metrics.csv")) |>
                       subset(algorithm == "RF", thr_value))   # max_sens_spec threshold of the RF (the only genus model)

# Continental Chile: Natural Earth 1:10m cropped to -76.5/-60, as in python/make_figures.py (Fig. 2)
world <- st_read(file.path(NE, "ne_10m_admin_0_countries.shp"), quiet = TRUE) |> st_make_valid()
chile_ll <- suppressWarnings(st_crop(world[world$ADMIN == "Chile", ], c(xmin = -76.5, xmax = -60, ymin = -60, ymax = -14)))
chile_a <- st_union(st_transform(chile_ll, albers))

# Rasters in Albers at 1 km (equal area per cell)
suit <- rast(file.path(R, "data/results/genus_main/ensemble_mean.tif"))
tpl  <- project(suit, albers, res = 1000)
suit <- mask(crop(tpl, vect(chile_a)), vect(chile_a)); names(suit) <- "suit"
mess <- project(rast(file.path(R, "mess/mess_genus.tif")), suit, method = "bilinear"); names(mess) <- "mess"
suitable_r <- (suit >= thr_rf) & (mess >= 0)

# Records (1,857 in continental Chile) and unique localities (coordinates rounded to 0.01 degrees)
reg <- read.csv("data/figshare/Syntrichia_occurrences_with_cleaning_flags.csv")   # 06_figshare_dataset.R
reg <- reg[toupper(as.character(reg$usedChileMaps)) == "TRUE", ]   # 1,857 records of continental Chile
stopifnot(nrow(reg) == 1857)
reg_sf <- st_as_sf(reg, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) |> st_transform(albers)
loc <- reg |> mutate(lon = round(decimalLongitude, 2), lat = round(decimalLatitude, 2)) |> distinct(lon, lat)
loc_sf <- st_as_sf(loc, coords = c("lon", "lat"), crs = 4326) |> st_transform(albers)
cat("Records:", nrow(reg), "| unique localities:", nrow(loc), "\n")

# ---------------------------------------------------------------- (1) hexagons
# Grid identical to that of Fig. 2 (fig_records_protected_areas.R): pointy-top hexagons,
# 50 km between opposite sides, origin at the limits of Chile 1:10m
hex_grid <- function(b, w) {
  R <- w / sqrt(3); dy <- 1.5 * R; polys <- list(); row <- 0; y <- b["ymin"] - R
  while (y < b["ymax"] + R) {
    x <- b["xmin"] - w + ifelse(row %% 2 == 1, w / 2, 0)
    while (x < b["xmax"] + w) {
      k <- 0:5; m <- cbind(x + R * sin(k * pi / 3), y + R * cos(k * pi / 3)); polys[[length(polys) + 1]] <- st_polygon(list(rbind(m, m[1, ])))
      x <- x + w
    }
    y <- y + dy; row <- row + 1
  }
  st_sf(geometry = st_sfc(polys, crs = albers))
}
hex <- hex_grid(st_bbox(chile_a), 50000)
hit <- unique(unlist(st_intersects(reg_sf, hex)))
hex <- hex[lengths(st_intersects(hex, chile_a)) > 0 | seq_len(nrow(hex)) %in% hit, ]
# Each record is assigned to a single hexagon (the first one), as in Fig. 2
assigned <- sapply(st_intersects(reg_sf, hex), `[`, 1)
hv <- vect(hex)
hex$mean_suitability <- terra::extract(suit, hv, fun = mean, na.rm = TRUE)[, 2]
hex$suitable_fraction <- terra::extract(suitable_r, hv, fun = mean, na.rm = TRUE)[, 2]
hex$analogous_fraction <- terra::extract(mess >= 0, hv, fun = mean, na.rm = TRUE)[, 2]
hex$n_records <- tabulate(assigned, nbins = nrow(hex))
hex$n_species <- sapply(seq_len(nrow(hex)), function(h) n_distinct(na.omit(reg$analysisSpecies[which(assigned == h)])))
cat('Hexagons:', nrow(hex), '| with records:', sum(hex$n_records > 0), '| max records:', max(hex$n_records), '| max richness:', max(hex$n_species), '| sum:', sum(hex$n_records), '\n')
hex$suitable <- !is.na(hex$suitable_fraction) & hex$suitable_fraction >= 0.5
hex$class <- with(hex, ifelse(!suitable, "Not suitable", ifelse(n_records == 0, "Suitable, no records",
                    ifelse(n_records <= 2, "Suitable, 1 or 2 records", "Suitable, 3+ records"))))
class_table <- hex |> st_drop_geometry() |> count(class) |> mutate(pct = round(100 * n / sum(n), 1))
print(class_table)
suit_hex <- hex[hex$suitable, ]
hex_summary <- data.frame(
  hex_total = nrow(hex), hex_suitable = nrow(suit_hex),
  suitable_without_records = sum(suit_hex$n_records == 0), pct_suitable_without_records = round(100 * mean(suit_hex$n_records == 0), 1),
  suitable_with_1_2 = sum(suit_hex$n_records %in% 1:2),
  pct_records_in_top10pct_hexagons = round(100 * sum(sort(hex$n_records, decreasing = TRUE)[1:ceiling(0.1 * nrow(hex))]) / sum(hex$n_records), 1),
  spearman_suitability_vs_records_in_suitable = round(cor(suit_hex$mean_suitability, suit_hex$n_records, method = "spearman"), 2))
print(hex_summary)
# Gaps by latitude band
cen <- st_coordinates(st_transform(st_centroid(suit_hex), 4326))
suit_hex$band <- cut(cen[, 2], c(-56, -45, -40, -35, -30, -25, -17), right = FALSE)
band_table <- suit_hex |> st_drop_geometry() |> group_by(band) |>
  summarise(hex_suitable = n(), without_records = sum(n_records == 0), pct_without = round(100 * mean(n_records == 0), 1),
            records = sum(n_records), .groups = "drop") |> arrange(desc(band))
print(band_table)

# ------------------------------------------------------- (2) accessibility bias
roads <- st_read(file.path(NE, "ne_10m_roads.shp"), quiet = TRUE) |> st_transform(4326)
roads <- suppressWarnings(st_crop(roads, c(xmin = -80, xmax = -62, ymin = -57, ymax = -16))) |> st_transform(albers)
towns <- st_read(file.path(NE, "ne_10m_populated_places.shp"), quiet = TRUE) |> st_transform(4326)
towns <- suppressWarnings(st_crop(towns, c(xmin = -80, xmax = -62, ymin = -57, ymax = -16))) |> st_transform(albers)
cat("Roads:", nrow(roads), "segments | populated places:", nrow(towns), "\n")
dist_km <- function(p, layer) as.numeric(st_distance(p, layer[st_nearest_feature(p, layer), ], by_element = TRUE)) / 1000

# Null models: (a) points with probability proportional to the suitability of the genus (all Chile)
#              (b) random points inside the suitable area (suitability >= threshold and analogous climate)
n_null <- 10000
dfi <- as.data.frame(suit, xy = TRUE, na.rm = TRUE)
null_a <- dfi[sample(nrow(dfi), n_null, replace = TRUE, prob = dfi$suit), c("x", "y")]
dfa <- as.data.frame(suitable_r, xy = TRUE, na.rm = TRUE); dfa <- dfa[dfa[, 3] == 1, ]
null_b <- dfa[sample(nrow(dfa), n_null, replace = TRUE), c("x", "y")]
as_sf <- function(d) st_as_sf(d, coords = c("x", "y"), crs = albers)

groups <- list(Localities = loc_sf, `Null suitability` = as_sf(null_a), `Null suitable area` = as_sf(null_b))
d <- bind_rows(lapply(names(groups), function(g) data.frame(group = g,
       road_km = dist_km(groups[[g]], roads), town_km = dist_km(groups[[g]], towns))))
accessibility_summary <- d |> group_by(group) |> summarise(n = n(),
  road_median_km = round(median(road_km), 1), pct_road_5km = round(100 * mean(road_km <= 5), 1),
  pct_road_10km = round(100 * mean(road_km <= 10), 1),
  town_median_km = round(median(town_km), 1), pct_town_25km = round(100 * mean(town_km <= 25), 1), .groups = "drop")
print(accessibility_summary)
tests <- bind_rows(lapply(c("Null suitability", "Null suitable area"), function(null) bind_rows(lapply(c("road_km", "town_km"), function(v) {
  x <- d[[v]][d$group == "Localities"]; y <- d[[v]][d$group == null]
  w <- wilcox.test(x, y, alternative = "less"); k <- suppressWarnings(ks.test(x, y))
  data.frame(null_model = null, variable = v, median_localities = round(median(x), 1), median_null = round(median(y), 1),
             median_ratio = round(median(x) / median(y), 2), wilcoxon_p = signif(w$p.value, 3),
             ks_D = round(unname(k$statistic), 3), ks_p = signif(k$p.value, 3))
}))))
print(tests)

# Save
write.csv(class_table, file.path(OUT, "hex_classes.csv"), row.names = FALSE)
write.csv(hex_summary, file.path(OUT, "hex_summary.csv"), row.names = FALSE)
write.csv(band_table, file.path(OUT, "hex_gaps_by_latitude.csv"), row.names = FALSE)
write.csv(accessibility_summary, file.path(OUT, "accessibility_summary.csv"), row.names = FALSE)
write.csv(tests, file.path(OUT, "accessibility_tests.csv"), row.names = FALSE)
saveRDS(list(hex = hex, d = d, roads = roads, towns = towns, loc = loc_sf, chile_a = chile_a, thr = thr_rf), file.path(OUT, "bias.rds"))
cat("RF threshold:", thr_rf, "\n")
