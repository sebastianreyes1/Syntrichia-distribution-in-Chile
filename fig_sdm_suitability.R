# ============================================================================
# fig_sdm_suitability.R
# Habitat suitability maps: genus Syntrichia + the five modelled species, Chile
# ----------------------------------------------------------------------------
# Figures produced (figures/):
#   Fig3_SDM_suitability.tiff / .pdf          Fig. 3: panels A-F, no records
#   FigS2_SDM_suitability_records.tiff / .pdf Fig. S2: same, calibration records overlaid
#   FigS3_MESS.tiff / .pdf                    Fig. S3: MESS surfaces      (MESS=1)
#   TableS_MESS_summary.csv                   area of novel climate        (MESS=1)
#   FigX_SDM_suitability_tgb.tiff / .pdf      target-group-background maps (TGB=1)
#   Fig_SDM_<taxon>.tiff                      one map per taxon (89 mm wide)
#
# Usage:
#   Rscript fig_sdm_suitability.R            # Fig. 3, Fig. S2 and single maps
#   MESS=1 Rscript fig_sdm_suitability.R     # Fig. S3 (needs 10_mess.R)
#   TGB=1  Rscript fig_sdm_suitability.R     # maps of the target-group-background models
#   SOLO_COMP=1 ...                          # skip the single-taxon maps
#
# Inputs:  data/results/<variant>_main/ensemble_mean.tif   (08_sdm_model.R)
#          data/clean/occ_<variant>.csv                      (05_clean_occurrences.R)
#          mess/mess_<model>.tif                             (10_mess.R)
# Requires: terra, sf, ggplot2, patchwork, rnaturalearthdata
#
# Style: Albers equal-area projection, same base map as Figs. 1, 2 and 4.
# Colour scale: viridis (low = dark purple, high = yellow), unsuitable = 0.
# ============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(patchwork)
})
sf_use_s2(FALSE)

use_tgb  <- Sys.getenv("TGB")  == "1"
use_mess <- Sys.getenv("MESS") == "1"

# Panel order = Fig. 3 (A genus, B-F species). names = file slug, values = results folder
taxa <- c(genus         = "genus_main",
          epilosa       = "epilosa_main",
          anderssonii   = "anderssonii_full_range_main",   # calibrated on its full range
          pseudorobusta = "pseudorobusta_main",
          breviseta     = "breviseta_main",
          squarripila   = "squarripila_main")
if (use_tgb) {   # sensitivity: species only
  taxa <- c(epilosa = "epilosa_tgb", anderssonii = "anderssonii_full_range_tgb",
            pseudorobusta = "pseudorobusta_tgb", breviseta = "breviseta_tgb",
            squarripila = "squarripila_tgb")
}
results_dir <- "data/results"; occ_dir <- "data/clean"; out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)

# ---- style ------------------------------------------------------------------
albers <- "+proj=aea +lat_1=-22 +lat_2=-50 +lat_0=-36 +lon_0=-71 +datum=WGS84 +units=m"
lon_lim <- c(-77, -64); lat_lim <- c(-56.3, -16.9)
col_sea <- "#bcd4e6"; col_neigh <- "#efe9df"; border_neigh <- "#d8cfc0"
col_chile <- "#ffffff"; border_chile <- "#6a6f73"; col_frame <- "#c9cdd1"
title_size <- 8.5; col_txt <- "#0b0b0b"; col_txt2 <- "#4a4f55"; col_axis <- "#7d8288"
pal <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")   # viridis
col_records <- "#ff5a5f"
font <- "Liberation Sans"   # on Windows use "Arial"
lw <- function(pt) pt / 2.845   # pt -> ggplot linewidth

# ---- base map ---------------------------------------------------------------
# rnaturalearthdata ships Natural Earth 1:50m. With rnaturalearthhires installed,
# use rnaturalearth::ne_countries(scale = 10, returnclass = "sf") for 1:10m.
countries <- rnaturalearthdata::countries50
if (!inherits(countries, "sf")) countries <- st_as_sf(countries)
box_ll <- st_as_sfc(st_bbox(c(xmin = -80, ymin = -58, xmax = -60, ymax = -15), crs = 4326))
countries <- suppressWarnings(st_intersection(st_make_valid(st_transform(countries, 4326)), box_ll))
chile <- countries[countries$admin == "Chile", ]
# Continental Chile only: drop Juan Fernandez and other oceanic islands (west of 77 W)
chile <- suppressWarnings(st_cast(chile, "POLYGON"))
chile <- chile[st_coordinates(st_centroid(chile))[, 1] > -77, ]
neighbours <- countries[countries$admin != "Chile", ]
chile_a <- st_transform(chile, albers); neigh_a <- st_transform(neighbours, albers)
chile_mask <- vect(st_union(chile_a))
frame_ll <- rbind(cbind(seq(lon_lim[1], lon_lim[2], length.out = 100), lat_lim[1]),
                  cbind(seq(lon_lim[1], lon_lim[2], length.out = 100), lat_lim[2]),
                  cbind(lon_lim[1], seq(lat_lim[1], lat_lim[2], length.out = 100)),
                  cbind(lon_lim[2], seq(lat_lim[1], lat_lim[2], length.out = 100)))
bb <- st_bbox(st_transform(st_sfc(st_multipoint(frame_ll), crs = 4326), albers))
x_lim <- c(bb["xmin"], bb["xmax"]); y_lim <- c(bb["ymin"], bb["ymax"])

pt_albers <- function(lon, lat) st_coordinates(st_transform(st_sfc(st_point(c(lon, lat)), crs = 4326), albers))
p_scale <- pt_albers(-70.3, -44.3)
p_scale_single <- pt_albers(-70.6, -40.6)   # single maps: legend sits at 45-48 S
p_north <- pt_albers(-65.2, -19.0)

# ---- rasters and records ------------------------------------------------------
read_suit <- function(folder) {
  r <- rast(file.path(results_dir, folder, "ensemble_mean.tif"))
  r <- project(r, albers, method = "bilinear", res = 1500)   # ~1.5 km, enough for the figure
  r <- mask(crop(r, chile_mask), chile_mask)
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(df)[3] <- "suitability"
  df
}
suit <- lapply(taxa, read_suit)

# Calibration records inside continental Chile (Argentine records are not shown)
records <- lapply(names(taxa), function(sp) {
  f <- file.path(occ_dir, paste0("occ_", if (sp == "anderssonii") "anderssonii_full_range" else sp, ".csv"))
  x <- read.csv(f)
  x <- st_transform(st_as_sf(x, coords = c("x", "y"), crs = 4326), albers)
  x[lengths(st_intersects(x, st_union(chile_a))) > 0, ]
})
names(records) <- names(taxa)

# ---- MESS (Elith et al. 2010), only if MESS=1 -----------------------------------
col_mess <- c("#dfe3e7", "#eb6834")
names(col_mess) <- c("Analogous climate (MESS ≥ 0)", "Novel climate (MESS < 0)")
mess_data <- NULL
if (use_mess) {
  mess_data <- lapply(names(taxa), function(sp) {
    tpl <- project(rast(file.path(results_dir, taxa[sp], "ensemble_mean.tif")), albers, res = 1500)
    r <- project(rast(sprintf("mess/mess_%s.tif", sp)), tpl, method = "bilinear")
    r <- mask(crop(r, chile_mask), chile_mask)
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(df)[3] <- "mess"
    df$class <- factor(ifelse(df$mess < 0, names(col_mess)[2], names(col_mess)[1]), levels = names(col_mess))
    df
  })
  names(mess_data) <- names(taxa)
}

# ---- one panel ----------------------------------------------------------------
panel <- function(sp, letter = NULL, with_records = FALSE, y_axis = TRUE, legend = TRUE,
                  scale_bar = TRUE, point_size = 1.1, scale_pos = p_scale,
                  legend_pos = c(0.80, 0.50), x_breaks = seq(-75, -65, 5), mess = FALSE) {
  ttl <- if (sp == "genus") {
    if (is.null(letter)) quote(italic("Syntrichia") ~ "(genus)" * phantom("p"))
    else bquote(bold(.(letter)) ~~ italic("Syntrichia") ~ "(genus)" * phantom("p"))
  } else {
    if (is.null(letter)) bquote(italic("S." ~ .(sp)) * phantom("p"))
    else bquote(bold(.(letter)) ~~ italic("S." ~ .(sp)) * phantom("p"))
  }
  g <- ggplot() +
    geom_sf(data = neigh_a, fill = col_neigh, colour = border_neigh, linewidth = lw(0.5)) +
    geom_sf(data = chile_a, fill = col_chile, colour = NA) +
    (if (mess) geom_raster(data = mess_data[[sp]], aes(x = x, y = y, fill = class))
     else      geom_raster(data = suit[[sp]],      aes(x = x, y = y, fill = suitability))) +
    geom_sf(data = chile_a, fill = NA, colour = border_chile, linewidth = lw(0.35)) +
    (if (mess) scale_fill_manual(values = col_mess, name = NULL, drop = FALSE)
     else scale_fill_gradientn(colours = pal, limits = c(0, 1), breaks = c(0, 0.5, 1),
                               labels = c("0", "0.5", "1"), name = "Habitat\nsuitability")) +
    labs(title = ttl, x = NULL, y = NULL)
  if (with_records)
    g <- g + geom_sf(data = records[[sp]], shape = 21, fill = col_records, colour = "white",
                     stroke = lw(0.25) * 2, size = point_size)
  if (scale_bar) {
    x0 <- scale_pos[1]; y0 <- scale_pos[2]
    g <- g + annotate("rect", xmin = x0, xmax = x0 + 250000, ymin = y0, ymax = y0 + 22000, fill = col_txt2) +
      annotate("rect", xmin = x0 + 250000, xmax = x0 + 500000, ymin = y0, ymax = y0 + 22000,
               fill = "white", colour = col_txt2, linewidth = lw(0.35)) +
      annotate("text", x = c(x0, x0 + 500000), y = y0 - 45000, label = c("0", "500 km"),
               size = 6 / .pt, family = font, colour = col_txt2)
  }
  g <- g + annotate("segment", x = p_north[1], xend = p_north[1], y = p_north[2] - 160000, yend = p_north[2],
                    arrow = arrow(length = unit(3.5, "pt"), type = "closed"),
                    linewidth = lw(0.6), colour = col_txt2) +
    annotate("text", x = p_north[1], y = p_north[2] + 55000, label = "N", size = 7 / .pt,
             family = font, fontface = "bold", colour = col_txt2)
  g + coord_sf(crs = albers, xlim = x_lim, ylim = y_lim, expand = FALSE, datum = 4326) +
    scale_x_continuous(breaks = x_breaks) + scale_y_continuous(breaks = seq(-55, -20, 5)) +
    theme_void(base_family = font, base_size = 7) +
    theme(panel.background = element_rect(fill = col_sea, colour = NA),
          panel.border = element_rect(fill = NA, colour = col_frame, linewidth = lw(0.5)),
          panel.grid = element_blank(),
          plot.title = element_text(size = title_size, colour = col_txt, hjust = 0, margin = margin(b = 1)),
          axis.text.y = if (y_axis) element_text(size = 6.7, colour = col_axis, hjust = 1, margin = margin(r = 2)) else element_blank(),
          axis.text.x = element_text(size = 6.7, colour = col_axis, margin = margin(t = 2)),
          legend.position = if (legend) legend_pos else "none",
          legend.background = element_rect(fill = "white", colour = col_frame, linewidth = lw(0.5)),
          legend.margin = margin(3, 4, 3, 4),
          legend.title = element_text(size = 6.5, colour = col_txt, lineheight = 0.9),
          legend.text = element_text(size = 6, colour = col_txt2),
          legend.key.width = unit(6, "pt"), legend.key.height = unit(22, "pt"),
          plot.margin = margin(3, 3, 3, 3))
}

save_fig <- function(g, name, width, height) {
  ggsave(file.path(out_dir, paste0(name, ".tiff")), g, width = width, height = height, dpi = 600,
         device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), g, width = width, height = height,
         device = cairo_pdf, bg = "white")
}
# In multi-panel figures the legend goes in a strip under the panels
bottom_legend <- theme(legend.position = "bottom", legend.direction = "horizontal",
                       legend.background = element_blank(),
                       legend.title = element_text(size = 9.5, colour = col_txt, vjust = 0.85, lineheight = 0.9),
                       legend.text = element_text(size = 9, colour = col_txt2),
                       legend.key.width = unit(28, "pt"), legend.key.height = unit(6, "pt"))

# ---- Fig. S3: MESS ---------------------------------------------------------------
if (use_mess) {
  summ <- do.call(rbind, lapply(names(taxa), function(sp) {
    m <- merge(mess_data[[sp]], suit[[sp]], by = c("x", "y"))
    high <- m$suitability >= 0.5
    data.frame(model = sp,
               pct_Chile_novel_climate = round(100 * mean(m$mess < 0), 1),
               km2_high_suitability = round(sum(high) * 2.25),              # 1.5 x 1.5 km cells
               pct_high_suitability_in_novel_climate = round(100 * mean(m$mess[high] < 0), 1),
               km2_high_suitability_analogous_climate = round(sum(high & m$mess >= 0) * 2.25))
  }))
  print(summ)
  write.csv(summ, file.path(out_dir, "TableS_MESS_summary.csv"), row.names = FALSE)
  ps <- lapply(seq_along(taxa), function(i)
    panel(names(taxa)[i], LETTERS[i], y_axis = i == 1, scale_bar = i == 1, point_size = 0.95,
          x_breaks = c(-75, -65), mess = TRUE))
  fig <- wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.background = element_blank(),
          legend.text = element_text(size = 9, colour = col_txt2),
          legend.key.width = unit(12, "pt"), legend.key.height = unit(9, "pt"))
  save_fig(fig, "FigS3_MESS", 7.2, 4.4)
  cat("Done (MESS). Figures in", out_dir, "\n"); quit(save = "no")
}

# ---- Fig. 3 (no records) and Fig. S2 (with records) ---------------------------------
for (with_rec in c(FALSE, TRUE)) {
  ps <- lapply(seq_along(taxa), function(i)
    panel(names(taxa)[i], LETTERS[i], with_records = with_rec, y_axis = i == 1, scale_bar = i == 1,
          point_size = 0.95, x_breaks = c(-75, -65)))
  fig <- wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect") & bottom_legend
  if (use_tgb) { if (!with_rec) save_fig(fig, "FigX_SDM_suitability_tgb", 7.2, 4.6); next }
  save_fig(fig, if (with_rec) "FigS2_SDM_suitability_records" else "Fig3_SDM_suitability", 7.2, 4.6)
}

# ---- one map per taxon (89 mm wide) ---------------------------------------------------
if (Sys.getenv("SOLO_COMP") != "1" && !use_tgb) {
  for (sp in names(taxa)) {
    g <- panel(sp, with_records = TRUE, point_size = 1.5, scale_pos = p_scale_single,
               legend_pos = c(0.80, 0.25))
    ggsave(file.path(out_dir, paste0("Fig_SDM_", sp, ".tiff")), g, width = 3.5, height = 7.4,
           dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  }
}
cat("Done. Figures in", out_dir, "\n")
