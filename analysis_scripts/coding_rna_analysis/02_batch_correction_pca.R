# FRDA Analysis Pipeline - 02 Batch Correction & PCA

# Load Config & Utils
source("analysis_scripts/coding_rna_analysis/00_init.R")
source("analysis_scripts/coding_rna_analysis/utils.R")

# Libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(limma)
    library(ggplot2)
    library(patchwork)
})

# Load Data
data_list <- readRDS(file.path(DIR_RDS, "data_clean.rds"))
counts <- data_list$counts
meta <- data_list$meta

# VST & Batch Correction
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = ~1)
vsd <- vst(dds, blind = FALSE)
expr_vst <- assay(vsd)

# Remove batch effects with respect to Condition, Sex, and Age
design_covs <- model.matrix(~ Condition + Sex + Age_scaled, data = meta)
expr_rb <- removeBatchEffect(expr_vst, batch = meta$Batch, design = design_covs)

# PCA Analysis
p_before <- plot_pca(expr_vst, meta,
    color_by = "Batch", shape_by = "Condition",
    title = "PCA Before Batch Correction"
)
p_after <- plot_pca(expr_rb, meta,
    color_by = "Batch", shape_by = "Condition",
    title = "PCA After Batch Correction"
)

# Save PCA plot
ggsave(file.path(DIR_PLOTS, "pca_batch_correction_all.png"),
    p_before + p_after,
    width = 16, height = 8, dpi = 600
)

# Save Corrected Data
saveRDS(
    list(vsd = vsd, expr_vst = expr_vst, expr_rb = expr_rb, meta = meta),
    file.path(DIR_RDS, "data_batch_corrected.rds")
)