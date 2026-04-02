# FRDA Analysis Pipeline - 08 Severity Trends
# Author: Senior Bioinformatician

# --- Load Config & Utils ---
source("analysis_scripts/00_init.R")
source("analysis_scripts/utils.R")

# --- Libraries ---
library(tidyverse)
library(limma)
library(splines)
library(edgeR)

# --- Load Data ---
meta_strat <- readRDS(file.path(DIR_RES, "meta_stratified.rds"))
data_corrected <- readRDS(file.path(DIR_RES, "data_corrected.rds"))
expr_vst <- data_corrected$expr_vst

# --- 1. Linear & Spline Regression ---
message(">>> Running Regression Analysis for Severity Traits...")
# Use patient-only counts for trajectory
counts_raw <- readRDS(file.path(DIR_RES, "data_clean.rds"))$counts
patients_idx <- which(meta_strat$Condition == "Patient")
meta_patients <- meta_strat[patients_idx, ]
counts_patients <- counts_raw[, patients_idx]

# Helper for Regression Loop
traits <- c("mFARS total", "FSA scale", "GAA1", "DD")
reg_results <- list()

for (tr in traits) {
    val <- as.numeric(as.character(meta_patients[[tr]]))
    if (all(is.na(val))) next

    # Filter patients with data
    keep <- !is.na(val)
    y_sub <- val[keep]
    counts_sub <- counts_patients[, keep]
    meta_sub <- meta_patients[keep, ]

    # Standardize trait
    y_scaled <- scale(y_sub)

    # Voom with Linear Model
    design_lin <- model.matrix(~y_scaled)
    dge <- DGEList(counts = counts_sub)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design_lin)
    fit <- lmFit(v, design_lin)
    fit <- eBayes(fit, robust = TRUE)
    res_lin <- topTable(fit, coef = 2, number = Inf) |> rownames_to_column("gene")

    # Voom with Spline Model (Non-linear)
    design_spl <- model.matrix(~ ns(y_scaled, df = 3))
    v_spl <- voom(dge, design_spl)
    fit_spl <- lmFit(v_spl, design_spl)
    fit_spl <- eBayes(fit_spl, robust = TRUE)
    res_spl <- topTable(fit_spl, coef = 2:4, number = Inf) |> rownames_to_column("gene")

    reg_results[[tr]] <- list(Linear = res_lin, Spline = res_spl)
}

# --- 2. Visualization of Trends ---
message(">>> Plotting Expression Trends...")
top_genes <- reg_results[["mFARS total"]]$Linear |>
    head(6) |>
    pull(gene)
p_trends <- plot_trends(meta_patients, expr_vst, top_genes, "mFARS total", "Top Genes vs mFARS Total")
ggsave(file.path(DIR_ML, "top_genes_mfars_trends.png"), p_trends, width = 12, height = 8)

# --- Save Results ---
saveRDS(reg_results, file.path(DIR_RES, "regression_severity_results.rds"))

message(">>> Severity Trends Analysis Complete.")
