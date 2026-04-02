suppressPackageStartupMessages({
    library(WGCNA)
    library(DESeq2)
    library(dplyr)
})

allowWGCNAThreads()

dds <- readRDS("ncrna_analysis/FRDA_ncRNA_filtered_dds.rds")
meta <- as.data.frame(colData(dds))

# VST
vsd <- vst(dds, blind = FALSE)
datExpr <- t(assay(vsd))

## Variance filter (top 10000)
vars <- apply(datExpr, 2, var)
top_genes <- names(sort(vars, decreasing = TRUE)[1:10000])
datExpr <- datExpr[, top_genes]

## Soft threshold
powers <- c(1:10, seq(12,20,2))

sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = "signed",
    verbose = 5
)

softPower <- sft$powerEstimate 
if (is.na(softPower) || softPower < 6) softPower <- 12
cat("Soft power:", softPower, "\n")

## Network construction
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

## Trait matrix
traits_df <- data.frame(
    Condition = as.numeric(meta$Condition == "Patient")
)

## Module–trait correlation
moduleColors <- labels2colors(net$colors)

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- bicor(MEs, traits_df, use = "p")
moduleTraitPvalue <- corPvalueStudent(
    moduleTraitCor,
    nSamples = nrow(datExpr)
)

## Save results
prefix <- paste0("ncrna_analysis/wgcna_patient_vs_control_ncrna")

write.csv(moduleTraitCor, paste0(prefix, "_cor.csv"))
write.csv(moduleTraitPvalue, paste0(prefix, "_pval.csv"))

gene_module <- data.frame(
    Gene = colnames(datExpr),
    Module = moduleColors
)

write.csv(gene_module,
            paste0(prefix, "_gene_modules.csv"),
            row.names = FALSE)

## Heatmap
png(paste0(prefix, "_heatmap.png"),
    width = 800, height = 1000)

labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = names(traits_df),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colors = blueWhiteRed(50),
    textMatrix = paste(
        signif(moduleTraitCor,2),
        "\n(",
        signif(moduleTraitPvalue,1),
        ")"
    ),
    zlim = c(-1,1),
    main = "Patient_vs_Control_ncRNA_WGCNA"
)

dev.off()
