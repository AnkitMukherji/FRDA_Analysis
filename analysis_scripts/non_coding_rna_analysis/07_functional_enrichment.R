suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(dplyr)
    library(ggplot2)
    library(readr)
})

# 1. LOAD DATA
val_targets <- readRDS("ncrna_analysis/Validated_miRNA_mRNA_Targets.rds")
# If RDS doesn't exist yet, we check for correlation CSVs
cor_all <- readRDS("ncrna_analysis/All_ncRNA_mRNA_Correlations.rds")

# 2. SEPARATE SETS FOR ENRICHMENT
# A: Validated miRNA Targets (Anti-correlated + Predicted)
mir_targets_mrs <- val_targets |>
    distinct(mRNA_ID) |>
    pull(mRNA_ID)

# B: lncRNA Co-expressed Hubs (Highly positive correlation > 0.6)
lnc_coexpressed_mrs <- cor_all |>
    filter(biotype != "miRNA", Correlation > 0.6) |>
    distinct(mRNA_ID) |>
    pull(mRNA_ID)

# 3. DEFINE THE ENRICHMENT FUNCTION
run_enrichment <- function(genes, label) {
    if (length(genes) < 10) {
        message("Too few genes for ", label, " (n = ", length(genes), "). Skipping.")
        return(NULL)
    }

    message("Running enrichment for ", label, " (n = ", length(genes), ")...")

    # Convert symbols to Entrez ID
    gene_ids <- bitr(genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
    )

    # GO Enrichment (BP - Biological Process)
    ego <- enrichGO(
        gene = gene_ids$ENTREZID,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
    )

    if (!is.null(ego) && nrow(ego) > 0) {
        p <- dotplot(ego, showCategory = 15) +
            ggtitle(paste("GO BP Enrichment:", label))
        ggsave(paste0("ncrna_analysis/GO_Enrichment_", label, ".png"), p, width = 8, height = 6)
        write.csv(as.data.frame(ego), paste0("ncrna_analysis/GO_Results_", label, ".csv"), row.names = FALSE)
    }

    # KEGG Enrichment
    ekegg <- enrichKEGG(
        gene = gene_ids$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 0.05
    )

    if (!is.null(ekegg) && nrow(ekegg) > 0) {
        p_k <- dotplot(ekegg, showCategory = 15) +
            ggtitle(paste("KEGG Enrichment:", label))
        ggsave(paste0("ncrna_analysis/KEGG_Enrichment_", label, ".png"), p_k, width = 8, height = 6)
        write.csv(as.data.frame(ekegg), paste0("ncrna_analysis/KEGG_Results_", label, ".csv"), row.names = FALSE)
    }

    return(list(GO = ego, KEGG = ekegg))
}

# 4. EXECUTE ENRICHMENT
# Enrichment for miRNA targets
tryCatch(
    {
        run_enrichment(mir_targets_mrs, "miRNA_Targets")
    },
    error = function(e) {
        message("Error in enrichment for miRNA targets: ", e$message)
    }
)

# Enrichment for lncRNA co-expression hubs
tryCatch(
    {
        run_enrichment(lnc_coexpressed_mrs, "lncRNA_Hubs")
    },
    error = function(e) {
        message("Error in enrichment for lncRNA hubs: ", e$message)
    }
)
