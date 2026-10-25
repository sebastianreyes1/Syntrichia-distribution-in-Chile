# =============================================================================
# 08_sdm_model.R  -  Syntrichia species distribution models (flexsdm)
#
# Usage (sbatch or console):
#   Rscript 08_sdm_model.R <variant> <mode>
#     variant: genus | epilosa | anderssonii | pseudorobusta | breviseta |
#              squarripila | anderssonii_full_range
#     mode:    test | main | tgb
#              test = quick test run
#              main = final models of the manuscript
#              tgb  = sensitivity analysis with target-group background
#                     (species only, not 'genus')
#
# Design of the models
#   - Input: data/clean/occ_<variant>.csv (output of 05_clean_occurrences.R)
#   - Species: calibration area = convex polygon of the presences with a 100 km
#     buffer (bmcp, Barve et al. 2011), inside the Chile + 200 km stack. All presences
#     fall inside the area, so all the folds are spatial. 10 pseudo-absences per
#     presence (maximum 10,000).
#   - Genus ('genus') and 'anderssonii_full_range': global background inside a 200 km
#     buffer around the presences (VanDerWal et al. 2009). The genus model uses as
#     many points as presences (prevalence 0.5), and the species full-range model
#     uses 10 times the presences (maximum 10,000). Spatial block folds over the whole
#     globe. Requires data/chelsa_global/CHELSA_bioXX.tif (script 02 or 07).
#   - Algorithms: GLM, GAM, MaxNet, RF and SVM (SVM is omitted for the genus model).
#   - Ensemble: only algorithms with Boyce > 0 and TSS >= 0.6.
#   - Final projection: Chile (stack_final_chile.tif).
#   - MaxNet: the final model is refitted with presences + pseudo-absences (section 4).
#     Without this, fit_max() called without 'background' leaves the final model
#     empty (0 coefficients). Cross-validation metrics are not affected.
#   - Mode 'tgb' (target-group background; Phillips et al. 2009) is a sensitivity
#     analysis. Everything else is identical to the main mode: same calibration
#     area, same number of points, same partition, same algorithms and same filter.
#     Only WHERE the background points come from changes: instead of random points,
#     cleaned records of other Syntrichia species (occ_genus.csv) are used, which
#     share the collection bias.
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args) >= 1) args[1] else "epilosa"
mode    <- if (length(args) >= 2) args[2] else "test"
stopifnot(mode %in% c("test", "main", "tgb"))
if (mode == "tgb" && variant == "genus") stop("Mode tgb is for species only")
# Variants calibrated with their full range (global background in a 200 km buffer
# around the presences). anderssonii_full_range: S. anderssonii with its whole
# subantarctic range (Chevalier et al. 2021, niche truncation).
global_range <- variant %in% c("genus", "anderssonii_full_range")

# ----------------------------- SETTINGS --------------------------------------
path_csv       <- file.path("data/clean", paste0("occ_", variant, ".csv"))
path_chile     <- "data/chile_layers/stack_final_chile.tif"
dir_global     <- "data/chelsa_global"
out_dir        <- file.path("data/results", paste0(variant, "_", mode))
vif_vars       <- c("bio02", "bio04", "bio08", "bio09", "bio15", "bio18", "bio19")
buffer_sp_m    <- 100000   # bmcp buffer, species
buffer_genus_m <- 200000   # background buffer, genus and full-range models
tss_threshold  <- 0.6
boyce_threshold <- 0
n_test         <- 2000     # presences of the genus model in test mode
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(flexsdm); library(terra); library(data.table); library(dplyr); library(sf)
})
set.seed(42)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out_dir, "log.txt")
msg <- function(...) { txt <- paste0(...); cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE) }

msg("==== Variant: ", variant, " | mode: ", mode, " | ", format(Sys.time()))

# Same scale transformation used when the climate was extracted at the occurrences
# (CHELSA temperatures / 10, see 02_extract_climate.R). Applying it to the global
# layers keeps presences, background and projection on the same scale.
rescale_chelsa <- function(r) {
  temp <- intersect(names(r), c("bio01","bio02","bio04","bio05","bio06","bio07",
                                "bio08","bio09","bio10","bio11"))
  for (v in temp) r[[v]] <- r[[v]] / 10
  r
}

# ==============================================================================
# 1. Cleaned presences
# ==============================================================================
pres <- as.data.frame(fread(path_csv))
pres <- pres[complete.cases(pres[, c("x", "y", vif_vars)]), c("x", "y", vif_vars)]
if (variant == "genus" && mode == "test" && nrow(pres) > n_test)
  pres <- pres[sample(nrow(pres), n_test), ]
pres$pr_ab <- 1L
msg("Presences: ", nrow(pres))

chile_stack <- rast(path_chile)
names(chile_stack) <- gsub(".*bio", "bio", names(chile_stack))
chile_stack <- chile_stack[[vif_vars]]

# Remove isolated fragments (Juan Fernandez / Rapa Nui)
patches <- terra::patches(chile_stack[[1]], directions = 8)
sizes   <- terra::freq(patches)
if (nrow(sizes) > 1)
  chile_stack <- terra::mask(chile_stack, patches == sizes$value[which.max(sizes$count)], maskvalue = FALSE)

# Target-group background: cleaned records of the genus (1 per 30" cell) inside the
# calibration area, excluding cells with presences of the species. If there are more
# than needed, a random sample is taken. If there are fewer, all are used and this
# is recorded in the log.
tgb_points <- function(env_M, n, env_layer = env_M) {
  tg <- as.data.frame(fread("data/clean/occ_genus.csv"))[, c("x", "y")]
  inside_M <- !is.na(terra::extract(env_M[[1]], as.matrix(tg))[, 1])
  tg <- tg[inside_M, ]
  tg_cell <- cellFromXY(env_layer, as.matrix(tg)); pr_cell <- cellFromXY(env_layer, as.matrix(pres[, c("x", "y")]))
  tg <- tg[!(tg_cell %in% pr_cell) & !duplicated(tg_cell), ]
  available <- nrow(tg)
  if (available > n) tg <- tg[sample(available, n), ]
  env <- terra::extract(env_layer, as.matrix(tg))[, vif_vars]
  out <- cbind(tg, env); out <- out[complete.cases(out), ]; out$pr_ab <- 0L
  msg("Target-group background: ", available, " genus records available in M | used: ", nrow(out),
      " (target ", n, ")")
  out
}

# ==============================================================================
# 2. Calibration area and pseudo-absences / background
# ==============================================================================
if (!global_range) {
  M <- calib_area(data = pres, x = "x", y = "y",
                  method = c("bmcp", width = buffer_sp_m), crs = crs(chile_stack))
  env_M <- crop(chile_stack, M, mask = TRUE)
  n_pa  <- min(10L * nrow(pres), 10000L)
  msg("Calibration area: bmcp + ", buffer_sp_m / 1000, " km | pseudo-absences: ", n_pa)
  if (mode == "tgb") {
    pa <- tgb_points(env_M, n_pa)
  } else {
    pa <- sample_pseudoabs(data = pres, x = "x", y = "y", n = n_pa, method = "random",
                           rlayer = env_M, calibarea = M)
    if (!all(vif_vars %in% names(pa))) {
      pa_xy <- pa[, c("x", "y")]; pa_xy$pr_ab <- 0L
      pa <- sdm_extract(data = pa_xy, x = "x", y = "y", env_layer = env_M, filter_na = TRUE)
    }
  }
  pa$pr_ab <- 0L
  env_part <- env_M            # layer for the spatial partition
} else {
  files <- file.path(dir_global, sprintf("CHELSA_bio%02d.tif", c(2, 4, 8, 9, 15, 18, 19)))
  if (!all(file.exists(files)))
    stop("Global layers are missing in ", dir_global, ". Run 07_download_chelsa_global.R")
  glob <- rast(files); names(glob) <- vif_vars
  glob <- rescale_chelsa(glob)

  # Aggregated version (~5', factor 10) used to sample the background and to
  # partition the data. The environmental values of the background are extracted
  # afterwards at 30".
  path_agg <- file.path(dir_global, "stack_global_agg10.tif")
  if (!file.exists(path_agg)) {
    msg("Aggregating global layers to ~5' (done only once)...")
    writeRaster(aggregate(glob, fact = 10, fun = "mean", na.rm = TRUE), path_agg, overwrite = TRUE)
  }
  glob_agg <- rast(path_agg)

  M <- calib_area(data = pres, x = "x", y = "y",
                  method = c("buffer", width = buffer_genus_m), crs = crs(glob_agg))
  env_M <- crop(glob_agg, M, mask = TRUE)
  # Genus: as many points as presences. Species: 10x presences (maximum 10,000),
  # as in the species-level models.
  n_bg  <- if (variant == "genus") nrow(pres) else min(10L * nrow(pres), 10000L)
  msg("Global background: ", buffer_genus_m / 1000, " km buffer around presences | n = ", n_bg)
  if (mode == "tgb") {
    pa <- tgb_points(env_M, n_bg, env_layer = glob)
  } else {
    bg_xy <- spatSample(env_M[[1]], size = n_bg, method = "random", na.rm = TRUE, xy = TRUE)[, c("x", "y")]
    bg_env <- terra::extract(glob, as.matrix(bg_xy))
    pa <- cbind(bg_xy, bg_env[, vif_vars]); pa$pr_ab <- 0L
  }
  env_part <- env_M
}

model_data <- bind_rows(pres[, c("x", "y", "pr_ab", vif_vars)],
                        pa[,   c("x", "y", "pr_ab", vif_vars)])
model_data <- model_data[complete.cases(model_data), ]
msg("Dataset: ", sum(model_data$pr_ab == 1), " presences / ", sum(model_data$pr_ab == 0), " absences")

# ==============================================================================
# 3. Spatial block partition (all points)
# ==============================================================================
inside <- complete.cases(terra::extract(env_part, model_data[, c("x", "y")], ID = FALSE))
if (any(!inside)) {
  msg("[WARNING] ", sum(!inside), " points have no value in the partition layer and are excluded")
  model_data <- model_data[inside, ]
}
data_copy <- unserialize(serialize(model_data, NULL))   # deep copy (works around a pr_ab bug in part_sblock)
part_out <- part_sblock(data = data_copy, env_layer = env_part, x = "x", y = "y", pr_ab = "pr_ab",
                        min_res_mult = 3L, max_res_mult = ifelse(global_range, 100L, 200L),
                        num_grids = ifelse(mode == "test", 10L, 30L), n_part = 5L, prop = 0.5)
part_col <- if (".part" %in% names(part_out$part)) ".part" else
            grep("part|fold|block", names(part_out$part), value = TRUE)[1]
if (nrow(part_out$part) != nrow(model_data)) stop("part_sblock returned a different number of rows")
model_data$fold <- as.integer(part_out$part[[part_col]])
msg("Presences per fold: ", paste(table(model_data$fold[model_data$pr_ab == 1]), collapse = " / "))
if (!"pr_ab" %in% names(model_data)) stop("model_data lost pr_ab during part_sblock")

# ==============================================================================
# 4. Individual models
# ==============================================================================
fit_model <- function(name, expr) {
  m <- tryCatch(eval(expr), error = function(e) { msg("  ERROR ", name, ": ", conditionMessage(e)); NULL })
  if (!is.null(m)) msg(sprintf("  %-7s TSS %.3f | AUC %.3f | Boyce %.3f", name,
                               m$performance$TSS_mean, m$performance$AUC_mean, m$performance$BOYCE_mean))
  m
}
args_fit <- list(data = model_data, response = "pr_ab", predictors = vif_vars,
                 partition = "fold", thr = "max_sens_spec")
models <- list(
  GLM    = fit_model("GLM",    quote(do.call(fit_glm, args_fit))),
  GAM    = fit_model("GAM",    quote(do.call(fit_gam, args_fit))),
  MaxNet = fit_model("MaxNet", quote(do.call(fit_max, c(args_fit, list(clamp = TRUE, pred_type = "cloglog"))))),
  RF     = fit_model("RF",     quote(do.call(fit_raf, args_fit))),
  # SVM is omitted in the genus model: the kernel matrix of kernlab::ksvm is O(n^2)
  # in memory and with ~40 thousand points it does not fit.
  SVM    = if (variant == "genus") NULL else fit_model("SVM", quote(do.call(fit_svm, args_fit))))
models <- Filter(Negate(is.null), models)

# Final MaxNet model: fit_max() evaluates each fold with presences and pseudo-absences,
# but fits the FINAL model dropping the absences and using 'background', which is not
# supplied here (the pseudo-absences go in 'data'). The final model is refitted with
# all the data and the same options as the folds ("default" classes, regmult = 1).
if (!is.null(models$MaxNet)) {
  f_mx <- maxnet::maxnet.formula(model_data$pr_ab, model_data[vif_vars], classes = "default")
  models$MaxNet$model <- maxnet::maxnet(p = model_data$pr_ab, data = model_data[vif_vars], f = f_mx,
                                        regmult = 1, addsamplestobackground = TRUE)
  msg("  Final MaxNet refitted with all the data: ", length(models$MaxNet$model$betas), " coefficients")
}
if (length(models) == 0) stop("No algorithm converged")

metrics <- bind_rows(lapply(names(models), function(nm)
  data.frame(algorithm = nm, variant = variant, models[[nm]]$performance)))
fwrite(metrics, file.path(out_dir, "model_metrics.csv"))
# Saved before the filter, so hours of computation are not lost if no model passes
# it (the decision can then be taken later without running again).
saveRDS(models,     file.path(out_dir, "individual_models.rds"))
saveRDS(model_data, file.path(out_dir, "model_dataset.rds"))

# ==============================================================================
# 5. Ensemble: Boyce > 0 and TSS >= 0.6
# ==============================================================================
passed <- metrics$algorithm[metrics$BOYCE_mean > boyce_threshold & metrics$TSS_mean >= tss_threshold]
if (mode == "tgb") {
  # Sensitivity: the SAME algorithms that entered the ensemble of the main run
  # (folder <variant>_main) are used, without filtering again. With target-group
  # background the metrics drop by construction (the background resembles the
  # presences, Phillips et al. 2009), so they are not comparable with those of the
  # main run. What is compared is the maps.
  met_main <- fread(file.path("data/results", paste0(variant, "_main"), "model_metrics.csv"))
  passed <- intersect(met_main$algorithm[met_main$BOYCE_mean > boyce_threshold & met_main$TSS_mean >= tss_threshold], names(models))
  msg("[tgb] Using the algorithms of the main ensemble: ", paste(passed, collapse = ", "))
}
msg("Models that pass the filter (Boyce > ", boyce_threshold, ", TSS >= ", tss_threshold, "): ",
    paste(passed, collapse = ", "))
if (length(passed) == 0) stop("No model passes the ensemble filter. Check before continuing.")
ensemble_models <- models[passed]

if (length(ensemble_models) >= 2) {
  m_ens <- fit_ensemble(models = ensemble_models, ens_method = "meanw",
                        thr_model = "max_sens_spec", metric = "TSS")
  fwrite(as.data.frame(m_ens$performance), file.path(out_dir, "ensemble_metrics.csv"))
} else {
  msg("[WARNING] Only one model passes the filter, it is used without an ensemble.")
  m_ens <- ensemble_models[[1]]
}

# ==============================================================================
# 6. Projection onto Chile
# ==============================================================================
flatten <- function(x) if (inherits(x, "SpatRaster")) x else rast(lapply(x, flatten))
preds <- rast(lapply(ensemble_models, function(m)
  flatten(sdm_predict(models = m, pred = chile_stack, thr = NULL, con_thr = FALSE))[[1]]))
writeRaster(mean(preds, na.rm = TRUE), file.path(out_dir, "ensemble_mean.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
writeRaster(stdev(preds, na.rm = TRUE), file.path(out_dir, "ensemble_sd.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
writeRaster(flatten(sdm_predict(models = m_ens, pred = chile_stack, thr = "max_sens_spec", con_thr = TRUE)),
            file.path(out_dir, "ensemble_flexsdm.tif"),
            overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))

saveRDS(models,     file.path(out_dir, "individual_models.rds"))
saveRDS(m_ens,      file.path(out_dir, "ensemble_model.rds"))
saveRDS(model_data, file.path(out_dir, "model_dataset.rds"))
msg("==== Done: ", out_dir, " | ", format(Sys.time()))
