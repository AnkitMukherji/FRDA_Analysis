suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(SummarizedExperiment)
  library(edgeR)
})

## 1. Load mRNA DESeq2 object and extract metadata
dds_mrna <- readRDS("dds_frda_all_ensg_removed.rds")

meta_final <- as.data.frame(colData(dds_mrna)) |>
  filter(Batch %in% c("8", "9")) |>
  rownames_to_column("sample_id")

## 2. Load ncRNA count matrix
ncrna_counts <- read.table(
  "ncrna_analysis/FRDA_Final_Annotated_Counts.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  quote = "",
  check.names = FALSE
)

## 3. Load RNAcentral mapping
id_mapping <- readRDS("ncrna_analysis/ids_mapping_rna_central_human.rds")
  
id_mapping_clean <- id_mapping |>
  filter(
    V4 == 9606,
    V5 %in% c("lncRNA","miRNA")
  ) |>
  mutate(
    priority = case_when(
      !is.na(V6) & V6 != "" & !str_detect(V6, "^ENS|^ENST|^URS|^HSAL") & V2 == "GENECARDS" ~ 1,
      !is.na(V6) & V6 != "" & !str_detect(V6, "^ENS|^ENST|^URS|^HSAL") & V2 == "HGNC" ~ 2,
      !is.na(V6) & V6 != "" & !str_detect(V6, "^ENS|^ENST|^URS|^HSAL") ~ 3,
      TRUE ~ 4
    )
  ) |>
  arrange(priority) |>
  transmute(
    rnacentral_id = V1,
    gene_symbol = V6,
    biotype = V5,
    source_db = V2,
    priority = priority
  )

lookup_urs <- id_mapping_clean |>
  distinct(rnacentral_id, .keep_all = TRUE)

lookup_symbol <- id_mapping_clean |>
  arrange(priority, rnacentral_id) |>
  group_by(gene_symbol) |>
  summarise(
    rnacentral_id = dplyr::first(rnacentral_id),
    biotype = dplyr::first(biotype),
    source_db = dplyr::first(source_db),
    .groups = "drop"
  )

## 4. Build feature annotation mapping
feature_tbl <- tibble(
  feature_id = rownames(ncrna_counts)
) |>
  mutate(
    id_type = if_else(str_detect(feature_id, "^URS"), "URS","SYMBOL")
  )

# 4a. URS mapping
feat_urs <- feature_tbl |>
  filter(id_type == "URS") |>
  left_join(lookup_urs,
            by = c("feature_id" = "rnacentral_id"))

# 4b. SYMBOL mapping
feat_symbol <- feature_tbl |>
  filter(id_type == "SYMBOL") |>
  left_join(lookup_symbol,
            by = c("feature_id" = "gene_symbol"))

# 4c. Merge both mappings and clean
feature_final <- bind_rows(feat_urs, feat_symbol) |>
  filter(!is.na(biotype)) |> 
  select(-gene_symbol, -priority) |>
  distinct(feature_id, .keep_all = TRUE) |>
  column_to_rownames("feature_id")

## 5. Subset counts and harmonize samples
genes_to_match <- intersect(rownames(ncrna_counts), rownames(feature_final))
counts_final <- ncrna_counts[genes_to_match, ]

# Rename sample columns (batch8.dragen.FA28.bam -> FA28)
colnames(counts_final) <- gsub(".*\\.([^.]+)\\.bam", "\\1", colnames(counts_final), ignore.case = TRUE)

# Match with metadata
shared_samples <- intersect(colnames(counts_final), meta_final$sample_id)
counts_final <- counts_final[, shared_samples]
meta_final <- meta_final[match(shared_samples, meta_final$sample_id), ]
rownames(meta_final) <- meta_final$sample_id

## 6. HIGH-CONFIDENCE EXPRESSION FILTERING
dge <- DGEList(counts = as.matrix(round(counts_final)))
keep <- filterByExpr(dge, group = meta_final$Condition)
counts_final <- counts_final[keep, ]
feature_final <- feature_final[rownames(counts_final), ]

message("High-confidence features retained: ", nrow(counts_final))

## 7. Save processed objects
dir.create("ncrna_analysis", showWarnings = FALSE)
saveRDS(feature_final, "ncrna_analysis/feature_ncrna.rds")
saveRDS(counts_final, "ncrna_analysis/counts_ncrna.rds")
saveRDS(meta_final, "ncrna_analysis/meta_ncrna.rds")