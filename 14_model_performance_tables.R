# 14_model_performance_tables.R
# Builds Table 1 (main text) and Table S2 (supplement) from the model outputs.
#
# Table S2: cross-validation performance of each algorithm in each model
#           (mean +/- SD across the 5 spatial blocks) and whether it was retained
#           in the ensemble (Boyce > 0 and TSS >= 0.6).
# Table 1:  records in Chile, presences used, algorithms retained and range of
#           TSS, AUC and Boyce across the retained algorithms of each ensemble.
#
# Inputs:  data/results/<variant>_main/model_metrics.csv and model_dataset.rds (08_sdm_model.R)
#          data/figshare/Chile_records_per_species.csv                          (06_figshare_dataset.R)
# Output:  figures/Table1_modelled_taxa.csv, figures/TableS2_model_performance.csv
suppressPackageStartupMessages(library(data.table))
dir.create("figures", showWarnings = FALSE)

tss_threshold <- 0.6; boyce_threshold <- 0   # same filter as 08_sdm_model.R
variants <- c(genus = "genus", epilosa = "epilosa", anderssonii = "anderssonii_full_range",
              pseudorobusta = "pseudorobusta", breviseta = "breviseta", squarripila = "squarripila")
fmt  <- function(m, s) sprintf("%.2f ± %.2f", m, s)
rng  <- function(x) if (length(x) == 0) "-" else if (round(min(x), 2) == round(max(x), 2)) sprintf("%.2f", x[1]) else sprintf("%.2f–%.2f", min(x), max(x))

s2 <- list(); t1 <- list()
for (nm in names(variants)) {
  dir_v <- file.path("data/results", paste0(variants[[nm]], "_main"))
  met <- fread(file.path(dir_v, "model_metrics.csv"))
  dat <- readRDS(file.path(dir_v, "model_dataset.rds"))
  n_pres <- sum(dat$pr_ab == 1); n_abs <- sum(dat$pr_ab == 0)
  met[, retained := BOYCE_mean > boyce_threshold & TSS_mean >= tss_threshold]
  s2[[nm]] <- met[, .(Model = if (nm == "genus") "Syntrichia (genus)" else paste0("S. ", nm),
                      Algorithm = algorithm, Presences = n_pres, `Pseudo-absences or background` = n_abs,
                      `TSS ± SD` = fmt(TSS_mean, TSS_sd), `AUC ± SD` = fmt(AUC_mean, AUC_sd),
                      `Boyce index ± SD` = fmt(BOYCE_mean, BOYCE_sd),
                      Retained = ifelse(retained, "Yes", "No"))]
  keep <- met[retained == TRUE]
  t1[[nm]] <- data.table(Taxon = if (nm == "genus") "Syntrichia (genus)" else paste0("S. ", nm),
                         `Presences used` = n_pres,
                         `Algorithms retained` = paste(keep$algorithm, collapse = ", "),
                         TSS = rng(keep$TSS_mean), AUC = rng(keep$AUC_mean), `Boyce index` = rng(keep$BOYCE_mean))
}
s2 <- rbindlist(s2); t1 <- rbindlist(t1)

# Records identified to each species in continental Chile (maps of Chile)
chile_sp <- fread("data/figshare/Chile_records_per_species.csv")
t1[, `Records in Chile` := ifelse(Taxon == "Syntrichia (genus)", sum(chile_sp$records),
                                  chile_sp$records[match(sub("^S\\. ", "Syntrichia ", Taxon), chile_sp$species)])]
setcolorder(t1, c("Taxon", "Records in Chile", "Presences used"))

fwrite(s2, "figures/TableS2_model_performance.csv")
fwrite(t1, "figures/Table1_modelled_taxa.csv")
print(t1); print(s2)
cat("Note: for the genus the presences are worldwide (thinned to one per 30\" cell);",
    "S. anderssonii was calibrated with its full range (outside Chile as well).\n")
