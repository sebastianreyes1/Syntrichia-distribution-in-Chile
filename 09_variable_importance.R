# =============================================================================
# 09_variable_importance.R  -  variable importance and response curves
# -----------------------------------------------------------------------------
# For each final model (folder <variant>_main) it uses ONLY the algorithms that
# passed the filter (TSS >= 0.6 and Boyce > 0) and computes:
#   - permutation importance (flexsdm::sdm_varimp): drop in TSS, AUC and Boyce
#     when each variable is permuted
#   - partial response curves (flexsdm::data_pdp): predicted suitability along
#     each variable, with the others at their mean
# Usage: Rscript 09_variable_importance.R <variant>
#        variant: genus | epilosa | anderssonii_full_range | pseudorobusta | breviseta | squarripila
# Output: data/results/<variant>_main/variable_importance.csv and response_curves.csv
# NOTE: temperatures (bio02, bio04, bio08, bio09) are on the scale of the model
#       (CHELSA / 10, see the README). They are converted back when plotting
#       (fig_variables.R), not here.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
variant <- if (length(args) >= 1) args[1] else "epilosa"
suppressPackageStartupMessages({ library(flexsdm); library(dplyr); library(data.table) })
set.seed(42)
vif_vars <- c("bio02", "bio04", "bio08", "bio09", "bio15", "bio18", "bio19")
dir_r  <- file.path("data/results", paste0(variant, "_main"))
models <- readRDS(file.path(dir_r, "individual_models.rds"))
model_data <- as.data.frame(readRDS(file.path(dir_r, "model_dataset.rds")))
met    <- fread(file.path(dir_r, "model_metrics.csv"))
passed <- met$algorithm[met$BOYCE_mean > 0 & met$TSS_mean >= 0.6]
models <- models[passed]
cat("Variant:", variant, "| algorithms:", paste(passed, collapse = ", "), "| rows:", nrow(model_data), "\n")
data_m <- model_data[, c("pr_ab", vif_vars)]

# ---- permutation importance ---------------------------------------------------------
n_sim <- if (variant == "genus") 5 else 20
imp <- bind_rows(lapply(names(models), function(nm) {
  r <- tryCatch(sdm_varimp(models = list(models[[nm]]), data = data_m, response = "pr_ab",
                           predictors = vif_vars, n_sim = n_sim, n_cores = 1,
                           thr = "max_sens_spec", clamp = TRUE, pred_type = "cloglog"),
                error = function(e) { cat("  ERROR varimp", nm, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) data.frame(algorithm = nm, as.data.frame(r)) else NULL
}))
fwrite(imp, file.path(dir_r, "variable_importance.csv"))
cat("Importance done:", nrow(imp), "rows\n")

# ---- partial response curves ----------------------------------------------------------
curves <- bind_rows(lapply(names(models), function(nm) bind_rows(lapply(vif_vars, function(v) {
  arg <- list(model = models[[nm]]$model, predictors = v, resolution = 100, training_data = data_m)
  if ("clamping" %in% names(formals(data_pdp))) arg$clamping <- TRUE   # depends on the flexsdm version
  r <- tryCatch(do.call(data_pdp, arg),
                error = function(e) { cat("  ERROR pdp", nm, v, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(NULL)
  d <- as.data.frame(r$pdpdata); names(d)[1:2] <- c("value", "suitability")
  data.frame(algorithm = nm, variable = v, d[, 1:2])
}))))
fwrite(curves, file.path(dir_r, "response_curves.csv"))
cat("Curves done:", nrow(curves), "rows\n==== Done\n")
