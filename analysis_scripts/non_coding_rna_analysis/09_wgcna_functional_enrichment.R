suppressPackageStartupMessages({
  library(WGCNA)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(tidyr)
})

# 1. LOAD DATA & METADATA
dds_nc <- readRDS("ncrna_analysis/FRDA_ncRNA_filtered_dds.rds")
dds_m <- readRDS("dds_frda_all_ensg_removed.rds")

# Shared samples
shared_samples <- intersect(colnames(dds_nc), colnames(dds_m))
dds_nc <- dds_nc[, shared_samples]
dds_m <- dds_m[, shared_samples]

# VST for both
vst_nc <- assay(vst(dds_nc, blind = FALSE))
vst_m <- assay(vst(dds_m, blind = TRUE))

# 2. LOAD WGCNA MODULES & SIGNIFICANCE
# Analysis A: Case-Control
res_a_pval <- read.csv("ncrna_analysis/wgcna_patient_vs_control_ncrna_pval.csv", row.names = 1)
res_a_mod <- read.csv("ncrna_analysis/wgcna_patient_vs_control_ncrna_gene_modules.csv")

# Analysis B: Patient-Only Clinical
# Check for modules significant in AT LEAST ONE primary trait
res_b_pval <- read.csv("ncrna_analysis/wgcna_patient_trait_pval.csv", row.names = 1)
res_b_mod <- read.csv("ncrna_analysis/wgcna_patient_gene_modules.csv")

# 3. IDENTIFY SIGNIFICANT MODULES (p < 0.05)
sig_mods_a <- rownames(res_a_pval)[res_a_pval$Condition < 0.05]
sig_mods_a <- gsub("^ME", "", sig_mods_a)

# For Analysis B, we check multiple clinical columns
clinical_traits <- c("mFARS_scaled", "mFARS_USS_scaled", "FSA_scaled")
sig_mods_b <- rownames(res_b_pval)[apply(res_b_pval[, clinical_traits], 1, function(x) any(x < 0.05))]
sig_mods_b <- gsub("^ME", "", sig_mods_b)

cat("Significant Case-Control modules:", length(sig_mods_a), "\n")
cat("Significant Patient-only Clinical modules:", length(sig_mods_b), "\n")

# 4. ENRICHMENT FUNCTION (Guilt-by-Association)
run_module_enrichment <- function(vst_m, vsd_nc, module_assignment, sig_mods, prefix) {
  dir.create(paste0("ncrna_analysis/enrichment_", prefix), showWarnings = FALSE, recursive = TRUE)
  
  results_list <- list()
  
  for (mod in sig_mods) {
    cat("Processing enrichment for module:", mod, "in", prefix, "\n")
    
    # Get ncRNAs in this module
    mod_nc <- module_assignment$Gene[module_assignment$Module == mod]
    if (length(mod_nc) == 0) next
    
    # Calculate Module Eigengene (ME) from expression
    # (Since we might not have the original ME matrix, we recalculate for precision)
    datExpr_mod <- t(vsd_nc[mod_nc, , drop = FALSE])
    me_mod <- moduleEigengenes(datExpr_mod, rep(mod, length(mod_nc)))$eigengenes
    
    # Correlate ME with all mRNAs (Bicor for robustness)
    cors <- bicor(me_mod[, 1], t(vst_m), use = "p")
    
    # Select top 500 positively correlated mRNAs (Guilt-by-Association)
    top_m <- names(sort(cors[1, ], decreasing = TRUE)[1:500])
    
    # Standardize gene names for clusterProfiler (assume ENSG or Symbols)
    # The dds_m has symbols as rownames but ENSG might be present
    genes_to_test <- top_m
    
    # GO Enrichment
    ego <- enrichGO(gene = genes_to_test,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'SYMBOL',
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)
    
    if (!is.null(ego) && nrow(ego) > 0) {
      write.csv(as.data.frame(ego), 
                paste0("ncrna_analysis/enrichment_", prefix, "/Module_", mod, "_GO_BP.csv"))
      
      p <- dotplot(ego, showCategory=15) + 
        ggtitle(paste0("Enrichment: Module ", mod, " (", prefix, ")"))
      ggsave(file.path("ncrna_analysis", paste0("enrichment_", prefix), paste0("Module_", mod, "_GO_dotplot.png")), 
             p, width = 8, height = 6)
      
      results_list[[mod]] <- ego
    }
  }
  return(results_list)
}

# 5. EXECUTE ENRICHMENT
# Run for Analysis A (Case-Control)
enrich_a <- run_module_enrichment(vst_m, vst_nc, res_a_mod, sig_mods_a, "CaseControl")

# Run for Analysis B (Patient Clinical)
# Note: Analysis B used only patient samples
vst_m_pat <- vst_m[, colnames(dds_nc)[dds_nc$Condition == "Patient"]]
vst_nc_pat <- vst_nc[, colnames(dds_nc)[dds_nc$Condition == "Patient"]]
enrich_b <- run_module_enrichment(vst_m_pat, vst_nc_pat, res_b_mod, sig_mods_b, "PatientClinical")

cat("Functional enrichment complete. Results saved in ncrna_analysis/enrichment_*/\n")
