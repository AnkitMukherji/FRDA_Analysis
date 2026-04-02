suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# 1. Load Prepared Data
counts_final <- readRDS("ncrna_analysis/counts_ncrna.rds")
meta_final <- readRDS("ncrna_analysis/meta_ncrna.rds")
feature_final <- readRDS("ncrna_analysis/feature_ncrna.rds")

# 2. Set Reference Level
meta_final$Condition <- relevel(meta_final$Condition, ref = "Control")

# 3. Construct DESeq2 Object
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_final)),
  colData = meta_final,
  design = ~ Batch + Age + Sex + Condition
)

dds <- DESeq(dds)

# 4. Differential Expression Results
res <- results(dds, contrast = c("Condition", "Patient", "Control"))
res_df <- as.data.frame(res) |>
  rownames_to_column("ncRNA_ID") |>
  left_join(
    feature_final |> rownames_to_column("ncRNA_ID"),
    by = "ncRNA_ID"
  ) |>
  arrange(padj)

# Save full results
write.csv(
  res_df,
  "ncrna_analysis/DESeq2_Patients_vs_Controls_All.csv",
  row.names = FALSE
)

feature_final <- feature_final[match(rownames(dds), rownames(feature_final)), ]
all(rownames(dds) == rownames(feature_final))
rowData(dds) <- cbind(
  rowData(dds),
  DataFrame(feature_final)
)
saveRDS(dds, "ncrna_analysis/FRDA_ncRNA_filtered_dds.rds")

# 5. Prepare Volcano Plot Data
vol_data <- res_df |>
  filter(!is.na(padj)) |>
  mutate(
    Expression = case_when(
      log2FoldChange >= 1  & padj <= 0.05 ~ "Upregulated",
      log2FoldChange <= -1 & padj <= 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    negLog10Padj = -log10(padj)
  )

# 6. Select Top Genes for Labeling
top_up <- vol_data |>
  filter(Expression == "Upregulated") |>
  filter(!grepl("ENSG|URS", ncRNA_ID)) |>
  arrange(padj) |>
  slice_head(n = 10)

top_down <- vol_data |>
  filter(Expression == "Downregulated") |>
  filter(!grepl("ENSG|URS", ncRNA_ID)) |>
  arrange(padj) |>
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)

# 7. Volcano Plot
p <- ggplot(vol_data,
            aes(x = log2FoldChange,
                y = negLog10Padj)) +
  geom_point(aes(color = Expression),
             alpha = 0.6,
             size = 1.8) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             color = "grey50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey50") +
  geom_text_repel(
    data = top_genes,
    aes(label = ncRNA_ID),
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.4
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "red3",
    "Downregulated" = "blue3",
    "Not significant" = "grey70"
  )) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Volcano Plot: Patients vs Controls (ncRNA)",
    x = "Log2 Fold Change",
    y = expression(-log[10]("Adjusted P-value")),
    color = ""
  )

# 8. Save Plot
ggsave(
  "ncrna_analysis/Volcano_Plot_Patients_vs_Controls.png",
  plot = p,
  width = 10,
  height = 5,
  dpi = 300
)