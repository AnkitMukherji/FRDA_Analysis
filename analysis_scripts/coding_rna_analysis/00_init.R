# FRDA Analysis Pipeline - Global Configuration

# Paths
DIR_RAW <- "salmon/"
DIR_DE_RES <- "deseq_de_results/"
DIR_VOLC_PLOTS <- "volcano_plot/"
DIR_PLOTS <- "plots/"
DIR_WGCNA <- "wgcna/"
DIR_PATHWAY <- "pathway_enrichment/"
DIR_ML <- "patient_stratification_prediction/"
DIR_RDS <- "rds_objects/"

# Ensure base directories exist
for (d in c(DIR_DE_RES, DIR_VOLC_PLOTS, DIR_WGCNA, DIR_PATHWAY, DIR_ML)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}