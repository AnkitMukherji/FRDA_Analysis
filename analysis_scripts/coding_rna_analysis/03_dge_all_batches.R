# FRDA Analysis Pipeline - 03 Differential Expression

# Load Config & Utils
source("analysis_scripts/coding_rna_analysis/00_init.R")
source("analysis_scripts/coding_rna_analysis/utils.R")

# Libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(tidyverse)
    library(writexl)
})

# Load Data
data_list <- readRDS(file.path(DIR_RDS, "data_clean.rds"))
counts <- data_list$counts
meta <- data_list$meta

expr_rb_list <- readRDS(file.path(DIR_RDS, "data_batch_corrected.rds"))
expr_rb <- expr_rb_list$expr_rb

# 1. Patient vs Control (All Batches)
design_all <- ~ Batch + Sex + Age_scaled + Condition
res_all <- dge(
    counts = counts, meta = meta, design = design_all,
    filename = file.path(DIR_DE_RES, "res_patient_vs_control_all.xlsx"),
    filename_volcano = file.path(DIR_VOLC_PLOTS, "volcano_patient_vs_control_all.png")
)

# 2. HCM Analysis
# Ensure HCM is 0/1 factor
meta_hcm <- meta
meta_hcm$HCM[is.na(meta_hcm$HCM)] <- 0
meta_hcm$HCM <- factor(meta_hcm$HCM, levels = c(0, 1))

# Patients Only for HCM comparison
meta_patients <- meta_hcm |> filter(Condition == "Patient")
counts_patients <- counts[, rownames(meta_patients)]

res_hcm_patients <- dge(
    counts = counts_patients, meta = meta_patients,
    design = ~ Batch + Sex + Age_scaled + HCM,
    filename = file.path(DIR_DE_RES, "res_hcm_vs_non_hcm_patients.xlsx"),
    filename_volcano = file.path(DIR_VOLC_PLOTS, "volcano_hcm_vs_non_hcm_patients.png")
)

# 3. FXN Expression Plot
fxn_data <- meta |>
    rownames_to_column("SampleID") |>
    mutate(FXN_Expr = expr_rb["FXN", ])

# Highlight specific samples from original script
highlight_samples <- c("FA18", "FA25", "FA43", "FA45")
fxn_data$Highlight <- ifelse(fxn_data$SampleID %in% highlight_samples, "Highlight", "Other")

p_fxn <- ggplot(fxn_data, aes(x = reorder(SampleID, FXN_Expr), y = FXN_Expr)) +
    geom_point(aes(color = Condition, alpha = Highlight, size = Highlight)) +
    scale_alpha_manual(values = c("Other" = 0.5, "Highlight" = 1)) +
    scale_size_manual(values = c("Other" = 2, "Highlight" = 4)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(title = "FXN Expression Across Samples", x = "Samples", y = "Batch-corrected VST Expression")

ggsave(file.path(DIR_PLOTS, "fxn_expression_plot.png"), p_fxn, width = 20, height = 8, dpi = 300)

# Save Results Summary
saveRDS(
    list(res_all = res_all, res_hcm_patients = res_hcm_patients),
    file.path(DIR_RDS, "dge_results_summary.rds")
)
