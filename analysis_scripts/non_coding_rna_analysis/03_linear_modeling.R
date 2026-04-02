suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(limma)
    library(ggplot2)
    library(dplyr)
    library(matrixStats)
})

dds <- readRDS("ncrna_analysis/FRDA_ncRNA_filtered_dds.rds")

# Subset patients only
dds_pat <- dds[, dds$Condition == "Patient"]
dds_pat$Condition <- droplevels(dds_pat$Condition)
meta <- as.data.frame(colData(dds_pat))

# Remove samples with zero total counts
keep_samps <- colSums(counts(dds_pat)) > 0
dds_pat <- dds_pat[, keep_samps]
meta <- meta[keep_samps, ]

# Feature Filtering (Top 10k variance)
dge_var <- DGEList(counts = counts(dds_pat))
dge_var <- calcNormFactors(dge_var)

log_cpm <- cpm(dge_var, log = TRUE)
vars <- rowVars(log_cpm)
names(vars) <- rownames(log_cpm)

top_10k <- names(sort(vars, decreasing = TRUE)[1:10000])

dge_subset <- dge_var[top_10k, ]
dge_subset <- calcNormFactors(dge_subset)

# Traits to analyze
traits <- c("FSA_scaled", "mFARS_scaled", "mFARS_USS_scaled")

for (trait in traits) {

    cat("\nRunning analysis for:", trait, "\n")

    ## metadata filtering
    valid_idx <- !is.na(meta[[trait]]) &
                 !is.na(meta$Age) &
                 !is.na(meta$Sex)

    meta_sub <- meta[valid_idx, ]
    dge_sub  <- dge_subset[, valid_idx]

    ## design matrix
    design <- model.matrix(
        as.formula(paste("~", trait, "+ Batch + Age + Sex")),
        data = meta_sub
    )

    ## voom + limma
    v <- voomWithQualityWeights(dge_sub, design, plot = FALSE)

    fit <- lmFit(v, design)
    fit <- eBayes(fit, robust = TRUE)

    ## results
    res <- topTable(
        fit,
        coef = trait,
        number = Inf,
        sort.by = "P"
    )

    res$GeneLabel <- rownames(res)

    cat("Significant (adj.P < 0.05):",
        sum(res$adj.P.Val < 0.05, na.rm = TRUE), "\n")

    cat("Significant (adj.P < 0.1):",
        sum(res$adj.P.Val < 0.1, na.rm = TRUE), "\n")

    ## save
    outfile <- paste0(
        "ncrna_analysis/FRDA_ncRNA_",
        trait,
        "_correlation.csv"
    )

    write.csv(res, outfile, row.names = FALSE)
}
