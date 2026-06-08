# Load Helper scripts
source("analysis_scripts/coding_rna_analysis/helper_scripts.R")

# Libraries
suppressPackageStartupMessages({
  library(readxl)
  library(DESeq2)
  library(tidyverse)
})

# Load RDS objects
rds_files <- list.files(path = "salmon", pattern = "\\.rds$", full.names = TRUE)
rds_list <- lapply(rds_files, readRDS)
names(rds_list) <- basename(rds_files)

# Metadata Preparation
attributes_frda <- read_xlsx("metadata/Metadata_RNAseq_FAvsHC.xlsx") |> subset(select = -Sno)
severity_scoring_frda <- read_xlsx("metadata/severity scoring.xlsx")
meta <- merge(attributes_frda, severity_scoring_frda, by = "Raw Data ID", all.x = TRUE) |>
    column_to_rownames("Raw Data ID")

meta[meta == "NA"] <- NA
meta <- prep_metadata(meta,
    factor_cols = c("Batch", "Sex", "Condition", "HCM", "Diabetes"),
    scale_cols = c("Age", "FSA scale", "mFARS total", "mFARS USS", "DD", "Onset", "GAA1", "GAA2")
)

# Standardize Batch column levels (remove .0)
meta$Batch <- factor(sub("\\.0$", "", as.character(meta$Batch)))

# Adding a separate column in metadata based on library prep method
meta <- meta |> mutate(Library_Prep = ifelse(as.numeric(as.character(Batch)) <= 7, "mRNA", "totalRNA"))

# Cleaning Samples
# Remove samples from batch 3 (rds_list[[3]]) that have been topped up in batch3_top (rds_list[[4]])
samples_batch3 <- colnames(assay(rds_list[[3]]))
samples_batch4 <- colnames(assay(rds_list[[4]]))
rds_list[[3]] <- rds_list[[3]][, !samples_batch3 %in% intersect(samples_batch3, samples_batch4)]

# Remove PC samples in Batch 8
rds_list[[9]] <- rds_list[[9]][, !colnames(assay(rds_list[[9]])) %in% c("PC1", "PC2", "PC3")]

# Combine Data
combined_rds <- do.call(cbind, rds_list)
counts <- assay(combined_rds)
counts <- counts[, sort(colnames(counts))]

# Map gene symbols to Ensembl IDs
gene_names <- rowData(rds_list[[1]])$gene_name
rownames(counts) <- make.unique(gene_names)

# Filter low-expression genes and remove genes with no annotation
counts_filt <- counts[(rowSums(counts >= 10) >= 5) & (!grepl("^ENSG", rownames(counts))), ]

# QC Selection for Batches 8 & 9
# Selecting samples that have million mapped reads to the transcriptome > 5 
# and million mapped fragments to the exonic region of the genome > 10
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
    file.path("rds_objects", "data_clean.rds")
)

# Also save for Batch 8/9 specific analysis
saveRDS(
    list(
        counts = counts_filt[, samples_qualified_89],
        meta = meta[samples_qualified_89, ]
    ),
    file.path("rds_objects", "data_batch89.rds")
)
