# Paper specific plots

# Sample-sample correlation of expression values across all genes
dds_sample_cor <- DESeqDataSetFromMatrix(countData = round(count_data_frda_all), 
                                          colData   = metadata_frda_all,
                                          design    = ~ 1) # no experimental variables are modeled or batch effects are removed
vst_mat_sample_cor <- vst(dds_sample_cor, blind = TRUE) # blind TRUE ignores design
cor_mat <- cor(assay(vst_mat_sample_cor), method = "pearson")
rb_mat_sample_cor <- removeBatchEffect(assay(vst_mat_sample_cor), batch = metadata_frda_all$Batch)
cor_mat_rb <- cor(rb_mat_sample_cor, method = "pearson")

batch_levels <- levels(metadata_frda_all$Batch)
batch_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(batch_levels)), batch_levels)
condition_colors <- c(Control = "#1F78B4", Patient = "#E31A1C")
ann_colors <- list(Condition = condition_colors, Batch = batch_colors)

annotation_col <- metadata_frda_all[, c("Condition", "Batch")]
annotation_col <- annotation_col[colnames(cor_mat), , drop = FALSE]

plot <- pheatmap(
  cor_mat_rb,
  name = "Pearson r",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "average",
  main = "Sample–Sample Correlation"
)
dev.copy(pdf, "QC/sample_sample_correlation_after_batch_effect_removal.pdf", width = 20, height = 20)
dev.off()

# Volcano plot with lfc < -1 and lfc > 1 with padj < 0.05 highlighted
res <- read_xlsx("res/res_all_frda_patients_vs_controls_ensg_removed.xlsx")
vol_data <- res |>
    na.omit() |> 
    mutate(Expression = case_when(
      log2FoldChange >= 1 & padj <= 0.05 ~ "Upregulated",
      log2FoldChange <= -1 & padj <= 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    ))

top_20 <- vol_data |>
  filter(Expression %in% c("Upregulated", "Downregulated")) |> 
  arrange(padj) |> 
  head(20)

p <- ggplot(vol_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Expression), alpha = 0.8, size = 1.5) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  geom_text_repel(
    data = top_20,
    aes(label = .data[["gene_names"]]),
    size = 1.5,
    max.overlaps = 20,
    force = 10, 
    segment.size = 0.2,
    segment.curvature = -0.5
  ) +
  scale_color_manual(values = c(
    "Upregulated" = "#CB477A",
    "Downregulated" = "#469FBF",
    "Not significant" = "gray70"
  )) +
  scale_x_continuous(name = expression(bold(log[2]("fold change")))) +
  scale_y_continuous(name = expression(bold(-log[10]("adj p")))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "black") +
  theme_classic(base_size = 3) +
  theme(
    axis.title = element_text(size = 6, face = "bold"),
    axis.text = element_text(size = 4),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    legend.position = c(0.8, 0.8),
    legend.justification = c(1, 1),
    legend.title = element_text(size = 5, face = "bold"),
    legend.text = element_text(size = 4),
    legend.key.size = unit(0.5, "line")
  )

ggsave("volcano_plot/volcano_all_frda_patients_vs_controls_final.png",
       plot = p,
       device = "png",
       units = "in",
       height = 8.5/3,
       width = 13/3,
       dpi = 600)

imp_pathways <- read_xlsx('pathway_enrichment/FRDA_All_Batches/imp_pathways.xlsx')
p <- imp_pathways |> 
  ggplot(aes(x = NES, y = reorder(ID, NES),
            fill = NES, alpha = -log10(p.adjust))) +
  geom_col(width = 0.7) +
  scale_fill_gradient2(
    low = "#3B8BC2",
    high = "#C93A5A",
    name = "NES"
  ) +
  scale_alpha_continuous(
    name = expression(-log[10]("adj p")),
    range = c(0.4, 1)) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  guides(
  fill = guide_colourbar(
    order = 1,
    barheight = unit(1.5, "cm"),
    barwidth  = unit(0.25, "cm"),
    title.theme = element_text(face = "bold", size = 5)
    ),
  alpha = guide_legend(
    order = 2,
    title = expression(bold(-log[10]("adj p"))),
    title.theme = element_text(face = "bold", size = 5),
    keywidth = unit(0.3, "cm"),
    keyheight = unit(0.3, "cm")
    )
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 6, face = "bold"),
    axis.text = element_text(size = 4, face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.3, colour = "grey40"),
    legend.position = "right",
    legend.text  = element_text(size = 4)
  )

ggsave("pathway_enrichment/FRDA_All_Batches/pathway_enrichment_final_test.png",
       plot = p,
       device = "png",
       units = "in",
       height = 8.5/3,
       width = 13/3,
       dpi = 600)

# GSVA and boxplot for immune cell types
immune_cell_genes <- read_csv("pathway_enrichment/FRDA_All_Batches/S5_gsva_list.csv")
immune_cell_genes_filt <- immune_cell_genes |> filter(symbol %in% rownames(vst_frda_all))
subset_vst <- vst_frda_all[immune_cell_genes_filt$symbol, ]
gene_sets_list <- split(
  immune_cell_genes_filt$symbol,
  immune_cell_genes_filt$cell_type
)

gsvapar <- gsvaParam(
      subset_vst,
      geneSets = gene_sets_list,
      assay = "counts",
      kcdf = "Gaussian",
      minSize = 1,
      maxDiff = TRUE,
      absRanking = FALSE
    )
es <- gsva(gsvapar)

immune_cells_boxplot_df <- assay(es) |>
  as.data.frame() |>
  rownames_to_column("Cell_type") |>
  pivot_longer(
    cols = -Cell_type,
    names_to = "Sample",
    values_to = "GSVA_score"
  ) |>
  left_join(
    metadata_frda_all |> rownames_to_column("Sample"),
    by = "Sample"
  )

stat_df <- immune_cells_boxplot_df |> 
  group_by(Cell_type) |> 
  wilcox_test(GSVA_score ~ Condition) |> 
  adjust_pvalue(method = "BH") |> 
  add_significance("p.adj") |> 
  left_join(
    immune_cells_boxplot_df |> 
      group_by(Cell_type) |> 
      summarise(y.position = max(GSVA_score, na.rm = TRUE) * 1.1),
    by = "Cell_type"
  )

immune_cells_boxplot_df_sig_nonsig <- immune_cells_boxplot_df |> 
  left_join(stat_df, by = "Cell_type")

plot <- ggplot(immune_cells_boxplot_df_sig_nonsig, aes(x = Condition, y = GSVA_score, fill = Condition)) +
  geom_boxplot(
    outlier.shape = NA, 
    alpha = 0.7, 
    width = 0.35,
    median.linewidth = 0.3
  ) +
  stat_pvalue_manual(
    stat_df,
    label = "p.adj.signif",
    size = 2,
    tip.length = 0.03,
    hide.ns = TRUE,
    y.position = "y.position"
  ) +
  scale_fill_manual(values = c(
      "Control" = "#3B8BC2",
      "Patient" = "#C93A5A"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  facet_wrap(~ Cell_type, scales = "free_y") +
  theme_pubclean(base_size = 7) +
  labs(
    title = NULL,
    x = NULL,
    y = "GSVA score"
    ) +
  theme(
    strip.text = element_text(size = 6, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 4),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 5, face = "bold"),
    legend.text = element_text(size = 4),
    legend.position = c(0.83, 0.18),
    # legend.box.background = element_rect(fill = "white", color = "black"),
    # legend.box.margin = margin(15, 20, 15, 20),  # top, right, bottom, left
    legend.background = element_blank())
ggsave("pathway_enrichment/FRDA_All_Batches/gsva_boxplot_immune_genes.png",
       plot = plot,
       device = "png",
       units = "in",
       height = 8.5/3,
       width = 13/3,
       dpi = 600)