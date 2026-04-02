# FRDA Analysis Pipeline - 01 Data Preprocessing

# Load Initial Config
source("analysis_scripts/coding_rna_analysis/00_init.R")
source("analysis_scripts/coding_rna_analysis/utils.R")

# Libraries
suppressPackageStartupMessages({
  library(readxl)
  library(DESeq2)
  library(tidyverse)
})

rds_files <- list.files(path = DIR_RAW, pattern = "\\.rds$", full.names = TRUE)
rds_list <- lapply(rds_files, readRDS)
names(rds_list) <- basename(rds_files)

# Metadata Preparation
attr_frda <- read_xlsx("Metadata_RNAseq_FAvsHC.xlsx") |> subset(select = -Sno)
sev_frda <- read_xlsx("severity scoring.xlsx")
meta <- merge(attr_frda, sev_frda, by = "Raw Data ID", all.x = TRUE) |>
    column_to_rownames("Raw Data ID")

meta[meta == "NA"] <- NA
meta <- prep_metadata(meta,
    factor_cols = c("Batch", "Sex", "Condition", "HCM", "Diabetes"),
    scale_cols = c("Age", "FSA scale", "mFARS total", "mFARS USS", "DD", "Onset", "GAA1", "GAA2")
)

# Library Prep Logic
# Standardize Batch column levels (remove .0)
meta$Batch <- factor(sub("\\.0$", "", as.character(meta$Batch)))

meta <- meta |>
    mutate(Library_Prep = ifelse(as.numeric(as.character(Batch)) <= 7, "mRNA", "totalRNA"))

# Cleaning Samples
# Remove batch3 top up samples
samples_batch3 <- colnames(assay(rds_list[[3]]))
samples_batch4 <- colnames(assay(rds_list[[4]]))
rds_list[[3]] <- rds_list[[3]][, !samples_batch3 %in% intersect(samples_batch3, samples_batch4)]

# Remove PCs in Batch 8
rds_list[[9]] <- rds_list[[9]][, !colnames(assay(rds_list[[9]])) %in% c("PC1", "PC2", "PC3")]

# Combine Data
combined_rds <- do.call(cbind, rds_list)
counts <- assay(combined_rds)
counts <- counts[, sort(colnames(counts))]

# Map Symbols
gene_names <- rowData(rds_list[[1]])$gene_name
rownames(counts) <- make.unique(gene_names)

# Filter low-expression and ENSG-only rows
counts_filt <- counts[(rowSums(counts >= 10) >= 5) & (!grepl("^ENSG", rownames(counts))), ]

# QC Selection for Batches 8 & 9
qc_result <- read_excel("QC/qc_dragen_salmon.xlsx", skip = 1)
samples_qualified_89 <- qc_result |>
    filter(Batch %in% c(8, 9)) |>
    mutate(Exonic = as.numeric(Exonic), M_Aligned = as.numeric(`M Aligned`)) |>
    filter(M_Aligned > 5 & Exonic > 10) |>
    pull(Sample)

# Final Selection
samples_b17 <- rownames(meta)[as.numeric(as.character(meta$Batch)) %in% 1:7]
samples_all <- c(samples_b17, samples_qualified_89)

meta_all <- meta[samples_all, ]
order_idx <- order(meta_all$Condition == "Patient", decreasing = TRUE)
meta_all <- meta_all[order_idx, ]
counts_all <- counts_filt[, samples_all][, order_idx]

# Save RDS
saveRDS(
    list(counts = counts_all, meta = meta_all, counts_raw = counts, meta_raw = meta),
    file.path(DIR_RDS, "data_clean.rds")
)

# Also save for Batch 8/9 specific analysis
saveRDS(
    list(
        counts = counts_filt[, samples_qualified_89],
        meta = meta[samples_qualified_89, ]
    ),
    file.path(DIR_RDS, "data_batch89.rds")
)
