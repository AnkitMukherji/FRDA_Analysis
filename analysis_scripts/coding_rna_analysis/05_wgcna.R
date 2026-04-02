# FRDA Analysis Pipeline - 05 WGCNA
# Author: Senior Bioinformatician

# --- Load Config & Utils ---
source("analysis_scripts/00_init.R")
source("analysis_scripts/utils.R")

# --- Libraries ---
library(WGCNA)
library(tidyverse)
library(DESeq2)
library(CorLevelPlot)
library(gridExtra)

# --- Load Data ---
data_corrected <- readRDS(file.path(DIR_RES, "data_corrected.rds"))
expr_rb <- data_corrected$expr_rb
meta <- data_corrected$meta

enableWGCNAThreads()

# --- 1. Gene Filtering for WGCNA ---
message(">>> Filtering genes for WGCNA...")
# Logic from .qmd line 1652: mean > 5 and top 50% variance
geneVars <- apply(expr_rb, 1, var)
keep_genes <- (rowMeans(expr_rb) > 5) & (geneVars > quantile(geneVars, 0.50))
expr_wgcna <- t(expr_rb[keep_genes, ])

# --- 2. Soft Thresholding ---
message(">>> Picking Soft Threshold...")
powers <- c(1:10, seq(12, 50, by = 2))
sft <- pickSoftThreshold(expr_wgcna, powerVector = powers, networkType = "signed", verbose = 3)

# Plot SFT
sft_data <- sft$fitIndices
p1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
    geom_text() +
    geom_hline(yintercept = 0.9, color = "red") +
    theme_classic()
p2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
    geom_text() +
    theme_classic()

png(file.path(DIR_WGCNA, "soft_threshold_selection.png"), width = 1000, height = 800, res = 150)
grid.arrange(p1, p2, nrow = 2)
dev.off()

# --- 3. Run WGCNA ---
message(">>> Building Co-expression Network...")
# Using soft power from init or sft
bwnet <- blockwiseModules(expr_wgcna,
    maxBlockSize = 20000, minModuleSize = 20,
    TOMType = "signed", power = WGCNA_POWER, mergeCutHeight = 0.15,
    randomSeed = 1234, verbose = 3
)

saveRDS(bwnet, file.path(DIR_WGCNA, "wgcna_result_all.rds"))

# --- 4. Module-Trait Correlation ---
message(">>> Calculating Module-Trait Correlations...")
# Prepare traits
meta_traits <- meta |>
    mutate(Condition_01 = ifelse(Condition == "Patient", 1, 0)) |>
    select(Condition_01, HCM, GAA1, GAA2, DD, `FSA scale`, `mFARS total`, `mFARS USS`) |>
    mutate(across(everything(), as.numeric))

# Fill NA
meta_traits[is.na(meta_traits)] <- 0

module_trait_cor <- cor(bwnet$MEs, meta_traits, use = "p")
module_trait_p <- corPvalueStudent(module_trait_cor, nrow(expr_wgcna))

# Heatmap
heatmap_data <- merge(bwnet$MEs, meta_traits, by = "row.names") |> column_to_rownames("Row.names")
png(file.path(DIR_WGCNA, "module_trait_heatmap.png"), width = 1800, height = 1200, res = 150)
CorLevelPlot(heatmap_data,
    x = colnames(meta_traits),
    y = colnames(bwnet$MEs),
    col = c("#398110ff", "#7ecc51ff", "white", "#f7afbeff", "#C93A5A")
)
dev.off()

# --- 5. Export Gene Lists ---
message(">>> Exporting gene lists for modules...")
gene_module_map <- data.frame(Gene = colnames(expr_wgcna), Module = bwnet$colors)
write_xlsx(gene_module_map, file.path(DIR_WGCNA, "gene_module_assignment.xlsx"))

message(">>> WGCNA Analysis Complete.")
