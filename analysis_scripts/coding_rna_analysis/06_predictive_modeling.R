# FRDA Analysis Pipeline - 06 Predictive Modeling
# Author: Senior Bioinformatician

# --- Load Config & Utils ---
source("analysis_scripts/00_init.R")
source("analysis_scripts/utils.R")

# --- Libraries ---
library(pROC)
library(caret)
library(tidyverse)
library(limma)

# --- Load Data ---
data_corrected <- readRDS(file.path(DIR_RES, "data_corrected.rds"))
vsd_obj <- data_corrected$vsd
meta <- data_corrected$meta

# Load DGE to find signature genes
res_all <- readRDS(file.path(DIR_RES, "dge_results_summary.rds"))$res_all

# --- 1. Stable Signature Genes Identification ---
message(">>> Finding stable signature genes...")
# Logic from .qmd line 2140: Synchronicity across batches and high significance
# For this script, we'll use top significant genes as a proxy signature
sig_genes <- res_all |>
    filter(padj < 0.01, abs(log2FoldChange) > 1.5) |>
    pull(gene_names)

expr_sig <- assay(vsd_obj)[sig_genes, ]
y <- meta$Condition

# --- 2. ROC Analysis for Individual Genes ---
message(">>> Running per-gene ROC...")
gene_roc <- calc_gene_roc_stats(expr_sig, y)
write_xlsx(gene_roc |> arrange(desc(AUC)), file.path(DIR_ML, "gene_level_ROC_results.xlsx"))

# --- 3. Cross-validated Centroid Scoring ---
message(">>> Running 10-fold CV for Centroid Score...")
set.seed(123)
folds <- createFolds(y, k = 10)
results_cv <- list()

for (k_size in c(20, 50, 100)) {
    cv_res <- map_dfr(seq_along(folds), function(i) {
        train_idx <- setdiff(seq_along(y), folds[[i]])
        test_idx <- folds[[i]]

        # Train
        expr_train <- expr_sig[, train_idx]
        y_train <- y[train_idx]

        # Select top k genes in training set
        gene_diff <- apply(expr_train, 1, function(g) median(g[y_train == "Patient"]) - median(g[y_train == "Control"]))
        top_k <- names(sort(abs(gene_diff), decreasing = TRUE))[1:min(k_size, length(gene_diff))]

        centroids <- train_centroids(expr_train[top_k, ], y_train)
        scores <- apply(expr_sig[top_k, test_idx], 2, signature_score, centroids = centroids)

        tibble(size = k_size, score = scores, truth = y[test_idx])
    })
    results_cv[[as.character(k_size)]] <- cv_res
}

# --- 4. Final Model & Threshold ---
message(">>> Building final predictive model...")
final_genes <- names(sort(abs(apply(expr_sig, 1, function(g) median(g[y == "Patient"]) - median(g[y == "Control"]))), decreasing = TRUE))[1:50]
final_centroids <- train_centroids(expr_sig[final_genes, ], y)
final_scores <- apply(expr_sig[final_genes, ], 2, signature_score, centroids = final_centroids)

roc_final <- roc(y, final_scores, levels = c("Control", "Patient"), direction = "auto")
threshold <- as.numeric(coords(roc_final, x = "best", best.method = "youden")["threshold"])

# Save Model
saveRDS(
    list(genes = final_genes, centroids = final_centroids, threshold = threshold),
    file.path(DIR_RES, "FRDA_diagnostic_signature_final.rds")
)

message(">>> Predictive Modeling Complete.")
