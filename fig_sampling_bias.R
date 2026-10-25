# Fig. S4: sampling gaps (hexagons) and accessibility bias (ECDF)
# Input:  figures/sampling_bias/bias.rds (12_sampling_bias.R)
# Output: figures/FigS4_sampling_bias.tiff and .pdf
suppressPackageStartupMessages({ library(sf); library(ggplot2); library(patchwork); library(dplyr) })
sf_use_s2(FALSE)
x <- readRDS("figures/sampling_bias/bias.rds"); hex <- x$hex; d <- x$d
albers <- "+proj=aea +lat_1=-22 +lat_2=-50 +lat_0=-36 +lon_0=-71 +datum=WGS84 +units=m"
font <- "Liberation Sans"; lw <- function(pt) pt / 2.845   # on Windows use "Arial"
col_sea <- "#bcd4e6"; col_nbr <- "#efe9df"; border_nbr <- "#d8cfc0"; border_chile <- "#6a6f73"
col_frame <- "#c9cdd1"; col_txt <- "#0b0b0b"; col_txt2 <- "#4a4f55"; col_axis <- "#7d8288"
countries <- st_as_sf(rnaturalearthdata::countries50) |> st_make_valid() |> st_transform(4326)
box_ll <- st_as_sfc(st_bbox(c(xmin = -80, ymin = -58, xmax = -60, ymax = -15), crs = 4326))
countries <- suppressWarnings(st_intersection(countries, box_ll))
neighbours_a <- st_transform(countries[countries$admin != "Chile", ], albers)
border_ll <- rbind(cbind(seq(-77, -64, length.out = 50), -56.3), cbind(seq(-77, -64, length.out = 50), -16.9),
                   cbind(-77, seq(-56.3, -16.9, length.out = 50)), cbind(-64, seq(-56.3, -16.9, length.out = 50)))
bb <- st_bbox(st_transform(st_sfc(st_multipoint(border_ll), crs = 4326), albers))
lev <- c("Suitable, 3+ records", "Suitable, 1 or 2 records", "Suitable, no records", "Not suitable")   # classes from 12_sampling_bias.R
cols <- c("#1c5cab", "#86b6ef", "#eb6834", "#ffffff")
hex$class <- factor(hex$class, levels = lev)
ma <- ggplot() +
  geom_sf(data = neighbours_a, fill = col_nbr, colour = border_nbr, linewidth = lw(0.5)) +
  geom_sf(data = x$chile_a, fill = "white", colour = NA) +
  geom_sf(data = hex, aes(fill = class), colour = "#d5dade", linewidth = lw(0.3)) +
  geom_sf(data = x$chile_a, fill = NA, colour = border_chile, linewidth = lw(0.35)) +
  scale_fill_manual(values = setNames(cols, lev), name = NULL) +
  guides(fill = guide_legend(ncol = 2)) +
  coord_sf(crs = albers, xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], expand = FALSE, datum = 4326) +
  scale_x_continuous(breaks = c(-75, -65)) + scale_y_continuous(breaks = seq(-55, -20, 5)) +
  labs(title = quote(bold("A") ~~ "Sampling gaps")) +
  theme_void(base_family = font, base_size = 7) +
  theme(panel.background = element_rect(fill = col_sea, colour = NA),
        panel.border = element_rect(fill = NA, colour = col_frame, linewidth = lw(0.5)),
        plot.title = element_text(size = 9, colour = col_txt, hjust = 0, margin = margin(b = 2)),
        axis.text.y = element_text(size = 6.7, colour = col_axis, hjust = 1, margin = margin(r = 2)),
        axis.text.x = element_text(size = 6.7, colour = col_axis, margin = margin(t = 2)),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white", colour = col_frame, linewidth = lw(0.5)),
        legend.margin = margin(3, 4, 3, 4), legend.text = element_text(size = 7, colour = col_txt2),
        legend.key.size = unit(8, "pt"))
d$group <- factor(d$group, levels = c("Localities", "Null suitability", "Null suitable area"),
                  labels = c("Collection localities (n = 885)", "Expected by suitability", "Random in suitable area"))
cg <- setNames(c("#eb6834", "#2a78d6", "#9aa3ab"), levels(d$group))
ecdf_plot <- function(v, title, xmax) ggplot(d, aes(x = .data[[v]], colour = group)) +
  stat_ecdf(linewidth = lw(1.2), pad = FALSE) +
  scale_colour_manual(values = cg, name = NULL) +
  coord_cartesian(xlim = c(0, xmax), expand = FALSE) +
  scale_y_continuous(labels = function(z) paste0(z * 100, "%"), breaks = seq(0, 1, 0.25)) +
  labs(title = title, x = "Distance (km)", y = "Cumulative proportion") +
  theme_minimal(base_family = font, base_size = 8) +
  theme(plot.title = element_text(size = 9, colour = col_txt, margin = margin(b = 2)),
        axis.text = element_text(colour = col_axis, size = 7), axis.title = element_text(colour = col_txt2, size = 7.5),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e1e0d9", linewidth = lw(0.4)),
        legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 7.5, colour = col_txt2),
        legend.key.height = unit(9, "pt"))
pb <- ecdf_plot("road_km", quote(bold("B") ~~ "Distance to the nearest main road"), 80)
pc <- ecdf_plot("town_km", quote(bold("C") ~~ "Distance to the nearest town"), 150)
fig <- ma + (pb / pc + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  plot_layout(widths = c(1, 1.35))
dir.create("figures", showWarnings = FALSE)
ggsave("figures/FigS4_sampling_bias.tiff", fig, width = 7.2, height = 7.0, dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave("figures/FigS4_sampling_bias.pdf", fig, width = 7.2, height = 7.0, device = cairo_pdf, bg = "white")
