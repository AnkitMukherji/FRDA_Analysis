# FRDA Analysis Pipeline - Utility Functions
# Description: Centralized helper functions for RNA-seq analysis, plots, and modeling.

# PREP
prep_metadata <- function(meta, factor_cols = NULL, scale_cols = NULL) {
    library(dplyr)
    m <- meta
    if (!is.null(factor_cols)) {
        for (c in factor_cols) {
            if (c %in% colnames(m)) {
                # Convert to character, remove .0 suffix if it exists, then factor
                # This ensures "1" and "1.0" become the same level
                val_char <- as.character(m[[c]])
                val_clean <- sub("\\.0$", "", val_char)
                m[[c]] <- factor(val_clean)
            }
        }
    }
    if (!is.null(scale_cols)) {
        for (c in scale_cols) {
            if (c %in% colnames(m)) {
                m[[paste0(c, "_scaled")]] <- scale(as.numeric(as.character(m[[c]])))
            }
        }
    }
    return(m)
}

# Check Design Matrix
# Check whether the model matrix is full rank or not 
# (i.e. there are no collinearities between covariates in the design formula)
check_design_collinearity <- function(formula, data) {
  design_matrix <- model.matrix(formula, data = data)
  # Check rank
  is_full_rank <- qr(design_matrix)$rank == ncol(design_matrix) # QR decomposition to check rank
  cat("Design matrix full rank: ", is_full_rank, "\n")
  if (!is_full_rank) {
    cat("\nDetected linear dependencies:\n")
    lm_obj <- lm(rep(1, nrow(design_matrix)) ~ design_matrix) # Fit linear model to identify dependencies
    print(alias(lm_obj)) # Display linear dependencies
  } else {
    cat("No linear dependencies detected.\n")
  }
  invisible(design_matrix)
}

# Plot: PCA
# PCA plot before and after removing batch effects with top 500 variable genes
# plotPCA(vsd, intgroup = c("Batch", "Condition")) # before batch effect removal
# or
plot_pca <- function(expr_matrix, metadata, top_n = 500, color_by, shape_by, title) {
    suppressPackageStartupMessages({
        library(ggplot2)
        library(matrixStats)
    })

    expr_matrix <- as.matrix(expr_matrix)
    metadata <- metadata[colnames(expr_matrix), , drop = FALSE]
    topVarGenes <- order(matrixStats::rowVars(expr_matrix), decreasing = TRUE)[seq_len(top_n)]
    
    pca_res <- prcomp(t(expr_matrix[topVarGenes, ]))
    percent_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
    
    pca_data <- data.frame(
        Sample = colnames(expr_matrix),
        PC1 = pca_res$x[, 1],
        PC2 = pca_res$x[, 2],
        Color = metadata[[color_by]]
    )
    
    if (!is.null(shape_by)) {
        pca_data$Shape <- metadata[[shape_by]]
    }
    
    pca_data <- na.omit(pca_data)
    
    p <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point(
        aes(color = Color, shape = if (!is.null(shape_by)) Shape else NULL),
        size = 3, alpha = 0.8
        ) +
        theme_bw(base_size = 14) +
        labs(
        title = title,
        x = paste0("PC1: ", round(percent_var[1] * 100, 1), "% variance"),
        y = paste0("PC2: ", round(percent_var[2] * 100, 1), "% variance"),
        color = color_by,
        shape = shape_by
        ) +
        theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
        )
    
    return(p)
}

# Plot: Volcano
volcano_plot <- function(res,
                         gene_label_col,
                         lfc_threshold = 2,
                         padj_threshold = 0.01,
                         top_n = 10,
                         base_size = 14) {
    
    suppressPackageStartupMessages({
        library(ggplot2)
        library(ggrepel)
        library(dplyr)
    })

    vol_data <- res |>
        na.omit() |> 
        mutate(Expression = case_when(
        log2FoldChange >= lfc_threshold & padj <= padj_threshold ~ "Upregulated",
        log2FoldChange <= -lfc_threshold & padj <= padj_threshold ~ "Downregulated",
        TRUE ~ "Not significant"
        ))

    top_up <- vol_data |>
        filter(Expression == "Upregulated") |>
        arrange(padj) |>
        head(top_n)
    
    top_down <- vol_data |>
        filter(Expression == "Downregulated") |>
        arrange(padj) |>
        head(top_n)
    
    top_genes <- bind_rows(top_up, top_down)

    p <- ggplot(vol_data, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(color = Expression)) +
        geom_text_repel(
        data = top_genes,
        aes(label = .data[[gene_label_col]]),
        size = 4,
        max.overlaps = 15
        ) +
        scale_color_manual(values = c(
        "Upregulated" = "#D55E00",
        "Downregulated" = "#0072B2",
        "Not significant" = "gray70"
        )) +
        scale_x_continuous(name = "Log2 fold change") +
        scale_y_continuous(name = "-Log10 adjusted p-value") +
        geom_hline(yintercept = -log10(padj_threshold), linetype = 2, color = "black") +
        geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = 2, color = "black") +
        theme_classic(base_size = base_size) +
        theme(
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        legend.position = "right"
        )
    
    return(p)
}

# DGE Wrapper
dge <- function(counts, meta, design, filename, filename_volcano,
                lfc_threshold = 2, padj_threshold = 0.01, top_n = 10) {
    
    suppressPackageStartupMessages({
        library(DESeq2)
        library(writexl)
        library(tibble)
    })

    dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                    colData = meta,
                                    design = design)
    dds <- DESeq(dds, parallel = FALSE)

    save_dds <- readline(prompt = "Save DESeqDataSet object? (yes/no): ")
    if (tolower(save_dds) == "yes") {
        dds_filename <- readline(prompt = "Enter filename to save DESeqDataSet object (e.g., dds.rds): ")
        saveRDS(dds, file.path(DIR_RDS, dds_filename))
    }
    
    cat("\nCoefficients:\n")
    print(resultsNames(dds))

    mode_choice <- readline(prompt = "\nUse 'coef' or 'contrast'? ")
    
    if (tolower(mode_choice) == "coef") {
        coef_choice <- readline(prompt = "Enter the coefficient name: ")
        res <- lfcShrink(dds, coef = coef_choice, type = "apeglm") |>
        as.data.frame() |>
        rownames_to_column("gene_names")
        
    } else if (tolower(mode_choice) == "contrast") {
        cat("\nEnter contrast as a single R expression.\n")
        cat("Examples:\n")
        cat("# Biological Question                               | Intercept Model (~ genotype + cond + g:c)                    | No-Intercept Model (~ 0 + group)\n")
        cat("# --------------------------------------------------------------------------------------------------------------------------\n")
        cat("# Condition effect in Genotype I (B vs A in I)       | results(dds, name=\"condition_B_vs_A\")                     | results(dds, contrast=list(\"groupI_B\", \"groupI_A\"))\n")
        cat("# Condition effect in Genotype II (B vs A in II)     | results(dds, list(c(\"condition_B_vs_A\", \"genotypeII.conditionB\"))) | results(dds, contrast=list(\"groupII_B\", \"groupII_A\"))\n")
        cat("# Interaction effect (Is B vs A different in II vs I?) | results(dds, name=\"genotypeII.conditionB\")               | results(dds, contrast=list(c(\"groupII_B\", \"groupI_A\"), c(\"groupII_A\", \"groupI_B\")))\n")
        cat("# ---------------------------------------------------------------------------\n")
        cat("Wrap the contrast in list(), e.g. list(\"Condition_Treated_vs_Control\") or list(\"Batch1.ConditionControl\", \"Batch9.ConditionControl\")\n\n")
        
        contrast_input <- readline(prompt = "Enter your contrast: ")
        contrast_list <- eval(parse(text = contrast_input))
        
        res <- lfcShrink(dds, contrast = contrast_list, type = "ashr") |>
        as.data.frame() |>
        rownames_to_column("gene_names")
        
    } else {
        stop("Invalid input")
    }

    resOrdered <- res[order(res$padj), ]
    write_xlsx(resOrdered, filename)

    volcano <- volcano_plot(resOrdered,
                            gene_label_col = "gene_names",
                            lfc_threshold = lfc_threshold,
                            padj_threshold = padj_threshold,
                            top_n = top_n)
    ggsave(
        filename = filename_volcano,
        plot = volcano,
        device = "jpeg",
        units = "in",
        height = 8,
        width = 15,
        dpi = 600
    )

    return(resOrdered)
}

# Analysis: GSEA
run_multi_gsea <- function(dge_df,
                           sets = c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:REACTOME",
                                    "TFT:GTRD", "TFT:TFT_LEGACY", "GO:BP", "GO:CC", "GO:MF", "IMMUNESIGDB"),
                           output_prefix = "pathway_enrichment",
                           species = "Homo sapiens") {
  
    suppressPackageStartupMessages({
        library(msigdbr)
        library(clusterProfiler)
        library(tidyverse)
        library(openxlsx)
        library(ggplot2)
        library(patchwork)
        library(enrichplot)
    })

    gene_sets <- lapply(setNames(sets, sets), function(x){
        msigdbr(species = species, subcollection = x)
    })
    
    # Filtering genes for those present in DESeq results
    gene_sets <- lapply(gene_sets, function(sets){
        sets[sets$gene_symbol %in% dge_df$gene_names, ]
    })

    filtered_df <- dge_df |> 
        filter(!is.na(log2FoldChange), !is.na(padj)) |>
        arrange(desc(abs(log2FoldChange))) |> 
        distinct(gene_names, .keep_all = TRUE)

    # Replacing padj = 0
    filtered_df <- filtered_df |> 
        mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj))

    # Ranking metric: sign(LFC) x –log10(padj)
    ranking_vector <- sign(filtered_df$log2FoldChange) * -log10(filtered_df$padj)
    names(ranking_vector) <- filtered_df$gene_names
    ranking_vector <- ranking_vector[is.finite(ranking_vector)]
    ranking_vector <- sort(ranking_vector, decreasing = TRUE)
    # Adding tiny jitter only if ties are present
    if (any(duplicated(ranking_vector))) {
        ranking_vector <- ranking_vector + rnorm(length(ranking_vector), 0, 1e-8)
        ranking_vector <- sort(ranking_vector, decreasing = TRUE)
    }

    gsea_results <- lapply(names(gene_sets), function(set_name){
        message(paste0("Running: ", set_name))
        
        gs_df <- gene_sets[[set_name]]
        
        gsea_res <- GSEA(
            geneList = ranking_vector,
            TERM2GENE = gs_df |> select(gs_name, gene_symbol),
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            eps = 0,
            nPermSimple = 10000,
            seed = TRUE
        )
        
        if (is.null(gsea_res) || nrow(gsea_res@result) == 0) {
        message(paste0("No significant pathways in ", set_name))
        return(NULL)
        }
        
        df <- gsea_res@result |> arrange(p.adjust)
        df$Description <- stringr::str_trunc(df$Description, width = 50)
        # Top pathway
        top_path <- df$ID[1]
        
        # Bubble plot
        if (nrow(df) >= 1) {
        p1 <- df |> 
            slice(1:min(20, nrow(df))) |>
            ggplot(aes(x = NES, y = reorder(Description, NES),
                    size = setSize, fill = p.adjust)) +
            geom_point(shape = 21) +
            theme_minimal() +
            scale_fill_viridis_c(trans = "log10") +
            ggtitle(paste0(set_name, " — Bubble Plot")) +
            xlab("NES") + ylab("")
        
        ggsave(
            paste0(output_prefix, "/GSEA_", gsub(":", "_", set_name), "_bubble.png"),
            p1, width = 10, height = 6
        )
        }
        
        # Running score plot
        p2 <- gseaplot2(
        gsea_res, geneSetID = top_path,
        title = paste0(set_name, ": ", df$Description[1])
        )
        
        ggsave(
        paste0(output_prefix, "/GSEA_", gsub(":", "_", set_name), "_running_score.png"),
        p2, width = 8, height = 6
        )
        
        return(gsea_res@result)
    })
    
    names(gsea_results) <- gsub(":", "_", sets)
    write.xlsx(gsea_results, paste0(output_prefix, "/GSEA_results.xlsx"))
    
    return(gsea_results)
}

# Analysis: GSVA
run_gsva_pipeline <- function(
  vst_object,
  msigdb_list = c("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:REACTOME",
                  "TFT:GTRD", "TFT:TFT_LEGACY", "GO:BP", "GO:CC", "GO:MF", "IMMUNESIGDB"),
  output_excel = "GSVA_results.xlsx",
  output_dir = "pathway_enrichment",
  assay_name = "counts",
  kcdf_mode = "Gaussian",
  width = 15,
  height = 10
) {

  suppressPackageStartupMessages({
    library(GSVA)
    library(SummarizedExperiment)
    library(msigdbr)
    library(limma)
    library(sva)
    library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    library(openxlsx)
    library(tidyverse)
    library(grid)
  })

  # ensure output dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  wb <- createWorkbook()

  message("Available columns in colData(vst_object):")
  print(colnames(colData(vst_object)))

  ## Ask user for design (variables used to build model)
  design_vars <- readline(prompt =
    "Enter variables for the design matrix (e.g. 'Batch + Sex + Age + Condition'): ")
  design_formula <- as.formula(paste0("~ ", design_vars))

  ## Ask user for annotation column (for heatmap)
  annot_col <- readline(prompt =
    "Enter the column name in colData(vst_object) to color in the heatmap (e.g. 'Condition'): ")
  if (!(annot_col %in% colnames(colData(vst_object)))) {
    stop("Annotation column not found in colData(vst_object). Aborting.")
  }
  annot_values <- unique(as.character(colData(vst_object)[[annot_col]]))
  message("Annotation groups detected:")
  print(annot_values)

  ## Ask user for annotation colors (produce named character vector)
  annot_colors <- character(length(annot_values))
  names(annot_colors) <- annot_values
  for (val in annot_values) {
    col_in <- readline(prompt = paste0("Enter color for group '", val,
                                       "' (e.g. 'salmon' or 'lightgreen') [press Enter to default 'grey']: "))
    if (col_in == "") col_in <- "grey"
    annot_colors[val] <- col_in
  }

  ## Build model matrices once from vst_object
  mod <- model.matrix(design_formula, colData(vst_object))
  mod0 <- model.matrix(~1, colData(vst_object))

  message("Available coefficients from the model (will be tested):")
  print(colnames(mod))
  coef_user <- readline(prompt = "Enter coefficient to test (exact column name, e.g. 'ConditionPatient'): ")
  if (!(coef_user %in% colnames(mod))) {
    hits <- grep(coef_user, colnames(mod), value = TRUE)
    if (length(hits) == 0) stop("Coefficient not found in model columns. Aborting.")
    message("Your input matched the following column(s): ", paste(hits, collapse = ", "))
    coef_user <- hits[1]
  }

  ## color scale for heatmap rows
  pwyexpcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  pwyexpcol <- rev(pwyexpcol)

  results_list <- list()

  for (collection in msigdb_list) {
    message("\nProcessing collection: ", collection, "\n")

    # cleaned name for files
    collection_clean <- gsub("[^A-Za-z0-9_]", "_", collection)
    collection_clean <- substr(collection_clean, 1, 31)

    msig_df <- tryCatch(
      msigdbr(species = "Homo sapiens", subcollection = collection),
      error = function(e1) tryCatch(msigdbr(species = "Homo sapiens", subcollection = collection),
                                   error = function(e2) NULL)
    )
    if (is.null(msig_df) || nrow(msig_df) == 0) {
      warning("msigdbr could not find collection '", collection, "'. Skipping.")
      next
    }

    gene_sets_list <- split(msig_df$gene_symbol, msig_df$gs_name)
    gene_sets_list <- lapply(gene_sets_list, function(pathway) intersect(pathway, rownames(vst_object)))
    gene_sets_list <- gene_sets_list[sapply(gene_sets_list, length) >= 1] # drop empty gene-sets

    if (length(gene_sets_list) == 0) {
      warning("No overlapping gene sets between collection '", collection, "' and vst_object. Skipping.")
      next
    }

    ## GSVA
    gsvapar <- gsvaParam(
      vst_object,
      geneSets = gene_sets_list,
      assay = assay_name,
      kcdf = kcdf_mode,
      minSize = 10,
      maxDiff = TRUE,
      absRanking = FALSE
    )

    es <- gsva(gsvapar)

    ## Save GSVA full matrix to workbook sheet
    gsva_df <- as.data.frame(assay(es)) |> rownames_to_column("Pathway")
    addWorksheet(wb, sheetName = collection_clean)
    writeData(wb, sheet = collection_clean, gsva_df)

    ## copy global model matrix and use locally
    mod_local <- mod
    mod0_local <- mod0

    ## SVA (estimate surrogate variables from GSVA scores) to estimate sample-level heterogeneity
    sv <- sva(assay(es), mod_local, mod0_local)
    mod_local <- cbind(mod_local, sv$sv)

    ## Fit limma
    fit <- lmFit(assay(es), mod_local)
    fit.eb <- eBayes(fit, robust = TRUE)
    gssizes <- geneSetSizes(es)
    fit.eb.trend <- eBayes(fit.eb, robust = TRUE, trend = gssizes)

    tt <- topTable(fit.eb.trend, coef = coef_user, n = Inf)
    DEpwys <- rownames(tt)[!is.na(tt$adj.P.Val) & tt$adj.P.Val <= 0.05]

    if (length(DEpwys) == 0) {
      message("No DE pathways found (FDR <= 0.05) for collection: ", collection)
      results_list[[collection_clean]] <- list(GSVA_matrix = assay(es), DE_table = tt, DE_matrix = NULL)
      next
    }

    message("Found ", length(DEpwys), " DE pathways for collection: ", collection)

    ## find condition columns in mod_local
    cond_cols <- which(colnames(mod_local) == coef_user)
    if (length(cond_cols) == 0) {
      ## try to infer variable name (e.g., ConditionPatient -> Condition)
      var_hint <- sub("([A-Za-z0-9_.-]+).*", "\\1", coef_user)
      cond_cols <- grep(paste0("^", var_hint), colnames(mod_local))
    }
    if (length(cond_cols) == 0) stop("Cannot determine condition column(s) for coef: ", coef_user)

    # covariates: drop intercept (col 1) and condition columns
    all_cols <- seq_len(ncol(mod_local))
    covar_cols <- setdiff(all_cols, c(1, cond_cols))
    if (length(covar_cols) == 0) covar_mat <- NULL else covar_mat <- mod_local[, covar_cols, drop = FALSE]

    # Remove covariates effect for plotting
    DEpwys_es_all <- removeBatchEffect(
      assay(es[DEpwys, , drop = FALSE]),
      covariates = covar_mat,
      design = mod_local[, cond_cols, drop = FALSE]
    )

    # Select top 30 pathways by adj.P.Val
    top_n <- min(30, nrow(DEpwys_es_all))
    top_idx <- order(tt[DEpwys, "adj.P.Val"])[1:top_n]
    DEpwys_es <- DEpwys_es_all[top_idx, , drop = FALSE]

    # Truncate pathway names
    rownames(DEpwys_es) <- stringr::str_trunc(rownames(DEpwys_es), width = 50)

    # cluster pathways
    gsetClust <- hclust(as.dist(1 - cor(t(DEpwys_es), method = "pearson")), method = "complete")

    # prepare annotation (use vst_object metadata)
    sample_annot <- as.character(colData(vst_object)[[annot_col]])
    names(sample_annot) <- colnames(es) # ensure names match columns in DEpwys_es

    ann <- HeatmapAnnotation(
      Annotation = sample_annot,
      col = list(Annotation = annot_colors),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )

    # truncate long pathway names for display (optional)
    rownames(DEpwys_es) <- stringr::str_trunc(rownames(DEpwys_es), width = 50)

    # Save heatmap
    out_pdf <- file.path(output_dir, paste0("GSVA_DEpws_", collection_clean, "_heatmap.pdf"))
    pdf(out_pdf, width = width, height = height)

    ht <- Heatmap(
      DEpwys_es,
      name = "GSVA Enrichment Score",
      col = pwyexpcol,
      cluster_rows = gsetClust,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 8),  
      column_names_gp = gpar(fontsize = 8),
      cluster_columns = FALSE,
      show_column_names = TRUE,
      top_annotation = ann,
      show_heatmap_legend = FALSE
    )

    gsva_enrichment_legend <- Legend(
      title = "GSVA Enrichment Score",
      col_fun = colorRamp2(
        breaks = c(-max(abs(DEpwys_es)), 0, max(abs(DEpwys_es))),
        colors = pwyexpcol[c(1, 128, 256)]
      ),
      at = seq(-max(abs(DEpwys_es)), max(abs(DEpwys_es)), length.out = 5),
      labels = round(seq(-max(abs(DEpwys_es)), max(abs(DEpwys_es)), length.out = 5), 2),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 5)
    )

    condition_legend <- Legend(
      labels = names(annot_colors),
      title = annot_col,
      legend_gp = gpar(fill = annot_colors),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )

    combine_legends <- packLegend(gsva_enrichment_legend, condition_legend, direction = "vertical")

    draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left", 
         heatmap_legend_list = combine_legends, padding = unit(c(5, 5, 5, 35), "mm")) #  bottom, left, top, right,
    dev.off()

    results_list[[collection_clean]] <- list(GSVA_matrix = assay(es), DE_table = tt, DE_matrix = DEpwys_es)

    message("Saved heatmap: ", out_pdf)
  }

  # save workbook inside output_dir
  out_wb <- file.path(output_dir, output_excel)
  saveWorkbook(wb, out_wb, overwrite = TRUE)
  message("Saved combined GSVA workbook: ", out_wb)

  message("GSVA pipeline finished. Returning results as a list.")
  return(invisible(results_list))
}