suppressPackageStartupMessages({
  library(WGCNA)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(SummarizedExperiment)
  library(matrixStats)
  library(tidyr)
})

allowWGCNAThreads()

# 1. LOAD DATA
counts_ncrna <- readRDS("ncrna_analysis/counts_ncrna.rds")
dds_mrna <- readRDS("dds_frda_all_ensg_removed.rds")
meta_shared <- readRDS("ncrna_analysis/meta_ncrna.rds")

# Match shared samples
shared_samples <- Reduce(
  intersect,
  list(
    colnames(dds_mrna),
    colnames(counts_ncrna),
    rownames(meta_shared)
  )
)

dds_shared <- dds_mrna[, shared_samples]
counts_ncrna <- counts_ncrna[, shared_samples]
meta_shared <- meta_shared[shared_samples, ]

# 2. VST TRANSFORMATION
vst_mrna <- assay(vst(dds_shared, blind = TRUE))

dds_ncrna <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_ncrna)),
  colData = meta_shared,
  design = ~1
)

vst_ncrna <- assay(vst(dds_ncrna, blind = TRUE))

# 3. SELECT SIGNIFICANT ncRNAs
de_ncrna <- read.csv("ncrna_analysis/DESeq2_Patients_vs_Controls_All.csv")

sig_nc <- de_ncrna |>
  filter(padj < 0.05) |>
  arrange(padj) |>
  pull(ncRNA_ID)

vst_ncrna_sig <- vst_ncrna[rownames(vst_ncrna) %in% sig_nc, , drop = FALSE]
message("Significant ncRNAs retained: ", nrow(vst_ncrna_sig))

# 4. REMOVE ZERO-VARIANCE GENES
keep_nc <- rowVars(vst_ncrna_sig) > 0
vst_ncrna_sig <- vst_ncrna_sig[keep_nc, ]

keep_m <- rowVars(vst_mrna) > 0
vst_mrna <- vst_mrna[keep_m, ]

# 5. VARIANCE FILTERING
var_cutoff <- quantile(rowVars(vst_mrna), 0.5)
vst_mrna <- vst_mrna[rowVars(vst_mrna) > var_cutoff, ]
message("mRNAs retained after filtering: ", nrow(vst_mrna))

# 6. PREPARE MATRICES FOR CORRELATION
df_nc <- t(vst_ncrna_sig) # samples × ncRNA
df_m <- t(vst_mrna) # samples × mRNA

# 7. CORRELATION ANALYSIS
cor_matrix <- cor(
  df_nc,
  df_m,
  method = "spearman",
  use = "everything"
)

# Compute p-values
p_matrix <- corPvalueStudent(cor_matrix, nrow(df_nc))

# 8. CONVERT TO LONG FORMAT
feature_info <- readRDS("ncrna_analysis/feature_ncrna.rds") |>
  rownames_to_column("ncRNA_ID")

cor_results <- as.data.frame(as.table(cor_matrix)) |>
  rename(
    ncRNA_ID = Var1,
    mRNA_ID = Var2,
    Correlation = Freq
  ) |>
  mutate(
    PValue = as.vector(p_matrix),
    FDR = p.adjust(PValue, method = "BH")
  ) |>
  left_join(feature_info, by = "ncRNA_ID")

# 9. FILTER SIGNIFICANT INTERACTIONS
# For miRNA, we prioritize anti-correlation (possible target repression)
# For lncRNA, we look at both (possible co-expression/guilt-by-association)
cor_results <- cor_results |>
  filter(FDR < 0.05) |>
  filter(
    (biotype == "miRNA" & Correlation < -0.3) |
      (biotype != "miRNA" & abs(Correlation) > 0.4)
  ) |>
  arrange(FDR)

message("Significant correlations: ", nrow(cor_results))

# 10. SAVE RESULTS
write.csv(
  head(cor_results, 5000),
  "ncrna_analysis/Top_ncRNA_mRNA_Correlations.csv",
  row.names = FALSE
)

saveRDS(
  cor_results,
  "ncrna_analysis/All_ncRNA_mRNA_Correlations.rds"
)
