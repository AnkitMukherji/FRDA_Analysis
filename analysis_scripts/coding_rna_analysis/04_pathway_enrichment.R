# FRDA Analysis Pipeline - 04 Pathway Enrichment

# Load Config & Utils
source("analysis_scripts/coding_rna_analysis/00_init.R")
source("analysis_scripts/coding_rna_analysis/utils.R")

# Libraries
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(tidyverse)
    library(GSVA)
    library(msigdbr)
    library(writexl)
})

# Load Data
res_all <- readRDS(file.path(DIR_RDS, "dge_results_summary.rds"))$res_all
data_corrected <- readRDS(file.path(DIR_RDS, "data_batch_corrected.rds"))
vst_obj <- data_corrected$vsd
meta <- data_corrected$meta

# 1. GO Enrichment (Overrepresentation analysis using hypergeometric test)
top_genes <- res_all |>
    filter(!is.na(padj) & (padj < 0.05) & (abs(log2FoldChange) >= 1)) |>
    pull(gene_names)

go_res <- enrichGO(
    gene = top_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL",
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05
)

write_xlsx(as.data.frame(go_res), file.path(DIR_PATHWAY, "FRDA_All_Batches/GO_Enrichment_All_Batches.xlsx"))

# Plot GO
p_go <- dotplot(go_res, showCategory = 15, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free_y")
ggsave(file.path(DIR_PATHWAY, "FRDA_All_Batches/GO_Dotplot_All_Batches.png"), p_go, width = 12, height = 15, dpi = 300)

# 2. GSEA (Multi-set)
gsea_frda_all <- run_multi_gsea(dge_df = dge_frda_all,
                                 sets = c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:REACTOME",
                                          "TFT:GTRD", "TFT:TFT_LEGACY", "GO:BP", "GO:CC", "GO:MF", "IMMUNESIGDB"),
                                 output_prefix = file.path(DIR_PATHWAY, "FRDA_All_Batches"),
                                 species = "Homo sapiens")

# 3. GSVA Pipeline
gsva_res <- run_gsva_pipeline(vst_object = vst_frda_all,
                                  msigdb_list = c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:REACTOME",
                                        "TFT:GTRD", "TFT:TFT_LEGACY", "GO:BP", "GO:CC", "GO:MF", "IMMUNESIGDB"),
                                  output_excel = "GSVA_results.xlsx",
                                  output_dir = file.path(DIR_PATHWAY, "FRDA_All_Batches"),
                                  assay_name = "counts",
                                  kcdf_mode = "Gaussian",
                                  width = 25,
                                  height = 10)