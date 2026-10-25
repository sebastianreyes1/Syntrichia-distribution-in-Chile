# =============================================================================
# 03_vif_selection.R
# Stepwise VIF selection (usdm::vifstep, threshold 5) of the 19 CHELSA variables
# at the occurrence records. The variables kept are bio02, bio04, bio08, bio09,
# bio15, bio18 and bio19.
#
# What a stepwise VIF does: for each candidate variable it fits a linear model
# against all the others and measures how much of its variance is explained by
# them (collinearity). If the highest VIF of the set is above the threshold
# (typically 5 or 10), that variable is removed and the calculation is repeated with
# the remaining ones, until all of them fall below the threshold.
#
# Input:  data/syntrichia_occ_with_climate.csv
# Output: data/syntrichia_occ_with_climate_vif.csv
# Requires: usdm, data.table
# =============================================================================
library(usdm)
library(data.table)

# --- 1. Load data and fix the duplicated-column problem ---
# Merging the sources left institutionCode_x / institutionCode_y (column name
# collision). We keep institutionCode_x (it comes from GBIF, the primary source),
# drop _y, and rename it to institutionCode.
dt <- fread("data/syntrichia_occ_with_climate.csv")

if ("institutionCode_x" %in% names(dt) && "institutionCode_y" %in% names(dt)) {
  dt[, institutionCode_y := NULL]
  setnames(dt, "institutionCode_x", "institutionCode")
}

# --- 2. Candidate variables ---
# The 19 CHELSA bioclimatic variables.
bio_cols <- sprintf("bio%02d", 1:19)
bio_cols <- bio_cols[bio_cols %in% names(dt)]

env <- dt[, ..bio_cols]
env <- na.omit(env)  # vifstep does not accept NA

cat("Rows used for VIF (no NA):", nrow(env), "of", nrow(dt), "\n")

# --- 3. Stepwise VIF ---
# VIF threshold = 5 (conservative for strongly correlated bioclimatic variables)
vif_threshold <- 5

vif_result <- vifstep(as.data.frame(env), th = vif_threshold)
print(vif_result)

selected_variables <- vif_result@results$Variables
cat("\nVariables kept (VIF <=", vif_threshold, "):\n")
print(selected_variables)

excluded_variables <- setdiff(bio_cols, selected_variables)
cat("\nVariables excluded because of collinearity:\n")
print(excluded_variables)

# --- 4. Save the final dataset with the selected variables only ---
final_cols <- c(setdiff(names(dt), bio_cols), selected_variables)
dt_final <- dt[, ..final_cols]

fwrite(dt_final, "data/syntrichia_occ_with_climate_vif.csv")
cat("\nDone: data/syntrichia_occ_with_climate_vif.csv (", nrow(dt_final), "rows,",
    ncol(dt_final), "columns ).\n")
