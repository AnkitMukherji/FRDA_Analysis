suppressPackageStartupMessages({
    library(WGCNA)
    library(DESeq2)
    library(dplyr)
    library(tibble)
})

allowWGCNAThreads()

# 1. LOAD DATA
dds <- readRDS("ncrna_analysis/FRDA_ncRNA_filtered_dds.rds")

# Subset patients only
dds_pat <- dds[, dds$Condition == "Patient"]
design(dds_pat) <- ~ 1
dds_pat$Condition <- droplevels(dds_pat$Condition)
meta_full <- as.data.frame(colData(dds_pat))

# 2. VST TRANSFORMATION
vsd <- vst(dds_pat, blind = FALSE)
datExpr <- t(assay(vsd))

# 3. SAMPLE AND FEATURE FILTERING
# Keep samples with at least one clinical trait
traits_of_interest <- c("mFARS_scaled", "mFARS_USS_scaled", "FSA_scaled")
valid_samps <- rowSums(!is.na(meta_full[, traits_of_interest])) > 0
datExpr <- datExpr[valid_samps, ]
meta <- meta_full[valid_samps, ]

cat("Samples used for consolidated WGCNA:", nrow(datExpr), "\n")

# Variance filter (top 10,000)
# Including more features after stringent expression filtering allows for capturing broader regulatory signals
vars <- apply(datExpr, 2, var)
top_genes <- names(sort(vars, decreasing = TRUE)[1:10000])
datExpr <- datExpr[, top_genes]

# 4. SOFT THRESHOLD SELECTION
powers <- c(1:10, seq(12,20,2))
sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = "signed",
    verbose = 5
)

softPower <- sft$powerEstimate 
if (is.na(softPower) || softPower < 6) softPower <- 12
cat("Selected Soft power:", softPower, "\n")

# 5. NETWORK CONSTRUCTION
net <- blockwiseModules(
    datExpr,
    power = softPower,
    corType = "bicor",
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3
)

# 6. TRAIT MATRIX
traits_df <- meta |>
    select(all_of(traits_of_interest), DD, GAA1) |>
    mutate(across(everything(), as.numeric)) |>
    mutate(Duration_x_GAA1 = DD * GAA1)

# 7. MODULE–TRAIT CORRELATION
moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Compute correlations for all traits at once
moduleTraitCor <- bicor(MEs, traits_df, use = "p")
moduleTraitPvalue <- corPvalueStudent(
    moduleTraitCor,
    nSamples = nrow(datExpr)
)

# 8. SAVE RESULTS
prefix <- "ncrna_analysis/wgcna_patient"

write.csv(moduleTraitCor, paste0(prefix, "_trait_cor.csv"))
write.csv(moduleTraitPvalue, paste0(prefix, "_trait_pval.csv"))

gene_module <- data.frame(
    Gene = colnames(datExpr),
    Module = moduleColors
)
write.csv(gene_module, paste0(prefix, "_gene_modules.csv"), row.names = FALSE)

# 9. HEATMAP
png(paste0(prefix, "_heatmap.png"), width = 1200, height = 1200, res = 150)

textMatrix <- paste0(
    signif(moduleTraitCor, 2),
    "\n(",
    signif(moduleTraitPvalue, 1),
    ")"
)
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = names(traits_df),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1,1),
    main = "Patient ncRNA WGCNA: Clinical Correlation"
)

dev.off()