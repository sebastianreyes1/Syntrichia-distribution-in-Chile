# 13_tgb_sensitivity.R
# Sensitivity to sampling bias: main map vs. target-group background, in continental
# Chile, overall and by zone (north < 30 S, centre-south >= 30 S).
# Inputs:  data/results/<model>_main/ensemble_mean.tif (main models, 08_sdm_model.R <variant> main)
#          data/results/<model>_tgb/ensemble_mean.tif  (08_sdm_model.R <variant> tgb)
# Output:  figures/TableS6_tgb_sensitivity.csv
suppressPackageStartupMessages({ library(terra); library(sf) })
sf_use_s2(FALSE)
dir.create("figures", showWarnings = FALSE)
p <- st_as_sf(rnaturalearthdata::countries50); p <- st_make_valid(p[p$admin == "Chile", ])
p <- suppressWarnings(st_cast(p, "POLYGON")); p <- p[st_coordinates(suppressWarnings(st_centroid(p)))[, 1] > -77, ]
chile <- vect(st_union(p))
sp <- c(epilosa = "epilosa", anderssonii = "anderssonii_full_range", pseudorobusta = "pseudorobusta",
        breviseta = "breviseta", squarripila = "squarripila")
D <- function(x, y) { px <- x / sum(x); py <- y / sum(y); 1 - 0.5 * sum(abs(px - py)) }   # Schoener's D
res <- do.call(rbind, lapply(names(sp), function(s) {
  a <- mask(crop(rast(sprintf("data/results/%s_main/ensemble_mean.tif", sp[s])), chile), chile)
  b <- mask(crop(rast(sprintf("data/results/%s_tgb/ensemble_mean.tif", sp[s])), chile), chile)
  xy <- crds(a, na.rm = FALSE); v <- cbind(values(a), values(b), xy[, 2]); v <- v[complete.cases(v), ]
  f <- function(i) c(r = cor(v[i, 1], v[i, 2], method = "spearman"), D = D(v[i, 1], v[i, 2]))
  t <- f(rep(TRUE, nrow(v))); n <- f(v[, 3] < -30)
  high <- v[, 1] >= 0.5
  data.frame(species = s, spearman_Chile = round(t["r"], 2), D_Chile = round(t["D"], 2),
             spearman_30_56S = round(n["r"], 2), D_30_56S = round(n["D"], 2),
             pct_high_main = round(100 * mean(high), 1), pct_high_tgb = round(100 * mean(v[, 2] >= 0.5), 1),
             pct_main_high_still_high_tgb = round(100 * mean(v[high, 2] >= 0.5), 1),
             tgb_suitability_mean_in_main_high = round(mean(v[high, 2]), 2), tgb_suitability_mean_elsewhere = round(mean(v[!high, 2]), 2))
}))
rownames(res) <- NULL; print(res); write.csv(res, "figures/TableS6_tgb_sensitivity.csv", row.names = FALSE)
