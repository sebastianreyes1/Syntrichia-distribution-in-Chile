# 10_mess.R
# MESS (Elith et al. 2010): similarity of the climate of every pixel of Chile to the
# calibration data (presences + background) of each model.
# Output: mess/mess_<model>.tif (MESS) and mess/mod_<model>.tif (most dissimilar variable)
#         <model> is genus, epilosa, anderssonii, pseudorobusta, breviseta or squarripila
#         (anderssonii uses its full-range model)
suppressPackageStartupMessages(library(terra))
dir.create("mess", showWarnings = FALSE)
v <- c("bio02","bio04","bio08","bio09","bio15","bio18","bio19")
s <- rast("data/chile_layers/stack_final_chile.tif"); names(s) <- gsub(".*bio", "bio", names(s)); s <- s[[v]]
models <- c(genus = "genus_main", epilosa = "epilosa_main", anderssonii = "anderssonii_full_range_main",
            pseudorobusta = "pseudorobusta_main", breviseta = "breviseta_main", squarripila = "squarripila_main")
cells <- which(complete.cases(values(s[[1]])))
X <- values(s)[cells, , drop = FALSE]
sim1 <- function(p, ref) {
  ref <- sort(ref); n <- length(ref); mn <- ref[1]; mx <- ref[n]
  f <- 100 * findInterval(p, ref, left.open = TRUE) / n          # % of ref < p
  out <- ifelse(f == 0, (p - mn) / (mx - mn) * 100,
         ifelse(f <= 50, 2 * f, ifelse(f < 100, 2 * (100 - f), (mx - p) / (mx - mn) * 100)))
  out
}
res <- list()
for (m in names(models)) {
  d <- readRDS(file.path("data/results", models[m], "model_dataset.rds"))
  S <- sapply(v, function(k) sim1(X[, k], d[[k]]))
  mess <- apply(S, 1, min); mod <- max.col(-S, ties.method = "first")
  r <- rast(s[[1]]); r[cells] <- mess; writeRaster(r, sprintf("mess/mess_%s.tif", m), overwrite = TRUE)
  r2 <- rast(s[[1]]); r2[cells] <- mod; writeRaster(r2, sprintf("mess/mod_%s.tif", m), overwrite = TRUE)
  res[[m]] <- data.frame(model = m, n_ref = nrow(d), pct_negative = round(100 * mean(mess < 0), 1),
                         limiting_variable = paste(names(sort(table(v[mod[mess < 0]]), decreasing = TRUE))[1:2], collapse = ", "))
  cat(m, "done\n")
}
print(do.call(rbind, res))
