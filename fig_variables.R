# Fig. S5: variable importance (permutation, drop in TSS) and response curves
# Inputs:  data/results/<model>_main/variable_importance.csv and response_curves.csv
#          (09_variable_importance.R, run for each model)
# Outputs: figures/FigS5_variable_importance_response.tiff and .pdf, figures/TableS4_variable_importance.csv
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
font <- "Liberation Sans"; lw <- function(pt) pt / 2.845   # on Windows use "Arial"
col_txt <- "#0b0b0b"; col_txt2 <- "#4a4f55"; col_axis <- "#7d8288"; grid_col <- "#e1e0d9"
mods <- c(genus = "Syntrichia (genus)", epilosa = "S. epilosa", anderssonii = "S. anderssonii",
          pseudorobusta = "S. pseudorobusta", breviseta = "S. breviseta", squarripila = "S. squarripila")
vars <- c(bio02 = "bio02\nDiurnal\nrange (°C)", bio04 = "bio04\nT seasonality\n(SD × 100)",
          bio08 = "bio08\nT wettest\nqtr (°C)", bio09 = "bio09\nT driest\nqtr (°C)",
          bio15 = "bio15\nP seasonality\n(CV, %)", bio18 = "bio18\nP warmest\nqtr (mm)",
          bio19 = "bio19\nP coldest\nqtr (mm)")
temp <- c("bio02", "bio04", "bio08", "bio09")   # on the model scale they are / 10 (see README): multiplied by 10 here

# Reads variable_importance.csv and response_curves.csv of each model (09_variable_importance.R)
folders <- c(genus = "genus_main", epilosa = "epilosa_main", anderssonii = "anderssonii_full_range_main",
             pseudorobusta = "pseudorobusta_main", breviseta = "breviseta_main", squarripila = "squarripila_main")
read_all <- function(file) do.call(rbind, lapply(names(folders), function(m)
  data.frame(read.csv(file.path("data/results", folders[m], file)), model = m)))
imp0 <- read_all("variable_importance.csv")
imp0 <- aggregate(cbind(TSS, AUC, BOYCE) ~ model + predictors, data = imp0, FUN = mean)
imp0$TSS_rel <- ave(pmax(imp0$TSS, 0), imp0$model, FUN = function(x) 100 * x / sum(x))
curves0 <- read_all("response_curves.csv")
dir.create("figures", showWarnings = FALSE)
imp <- imp0 |>
  mutate(model = factor(model, levels = names(mods), labels = mods),
         variable = factor(predictors, levels = names(vars)))
ha <- ggplot(imp, aes(x = variable, y = model, fill = TSS_rel)) +
  geom_tile(colour = "white", linewidth = lw(1)) +
  geom_text(aes(label = sprintf("%.0f", TSS_rel), colour = TSS_rel > 18), size = 7 / .pt, family = font) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = col_txt), guide = "none") +
  scale_fill_gradientn(colours = c("#f4f8fd", "#cde2fb", "#86b6ef", "#2a78d6", "#0d366b"), limits = c(0, 35),
                       name = "Relative\nimportance (%)") +
  scale_y_discrete(limits = rev(mods)) +
  scale_x_discrete(labels = names(vars)) +
  labs(title = quote(bold("A") ~~ "Variable importance (drop in TSS after permutation)"), x = NULL, y = NULL) +
  theme_minimal(base_family = font, base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "italic", colour = col_txt2, size = 7.5),
        axis.text.x = element_text(colour = col_txt2, size = 7.5), plot.title = element_text(size = 9, colour = col_txt),
        legend.title = element_text(size = 7.5, colour = col_txt), legend.text = element_text(size = 7, colour = col_txt2),
        legend.key.height = unit(10, "pt"), legend.key.width = unit(6, "pt"))

curves <- curves0 |>
  mutate(value = ifelse(variable %in% temp, value * 10, value)) |>
  group_by(model, variable, x = round(value, 6)) |>
  summarise(mean_suit = mean(suitability), mn = min(suitability), mx = max(suitability), .groups = "drop") |>
  mutate(model = factor(model, levels = names(mods), labels = mods),
         variable = factor(variable, levels = names(vars), labels = vars))
cb <- ggplot(curves, aes(x = x)) +
  geom_ribbon(aes(ymin = mn, ymax = mx), fill = "#cde2fb") +
  geom_line(aes(y = mean_suit), colour = "#1c5cab", linewidth = lw(1)) +
  facet_grid(model ~ variable, scales = "free_x", switch = "x") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(n.breaks = 3) +
  labs(title = quote(bold("B") ~~ "Partial response curves (mean and range of the retained algorithms)"),
       x = NULL, y = "Predicted suitability") +
  theme_minimal(base_family = font, base_size = 7) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid_col, linewidth = lw(0.4)),
        strip.text.y = element_text(face = "italic", angle = 0, hjust = 0, colour = col_txt2, size = 7),
        strip.text.x = element_text(colour = col_txt2, size = 6.5, lineheight = 0.9), strip.placement = "outside",
        axis.text = element_text(colour = col_axis, size = 5.5), axis.title.y = element_text(colour = col_txt2, size = 7.5),
        panel.spacing.x = unit(4, "pt"), panel.spacing.y = unit(4, "pt"),
        plot.title = element_text(size = 9, colour = col_txt))
fig <- ha / cb + plot_layout(heights = c(1.25, 5))
ggsave("figures/FigS5_variable_importance_response.tiff", fig, width = 7.2, height = 9.0, dpi = 600, device = "tiff", compression = "lzw", bg = "white")
ggsave("figures/FigS5_variable_importance_response.pdf", fig, width = 7.2, height = 9.0, device = cairo_pdf, bg = "white")
# Table S4
tab <- reshape(imp[, c("model", "variable", "TSS_rel")], idvar = "model", timevar = "variable", direction = "wide")
names(tab) <- sub("TSS_rel.", "", names(tab)); tab[, -1] <- round(tab[, -1], 1)
write.csv(tab, "figures/TableS4_variable_importance.csv", row.names = FALSE)
