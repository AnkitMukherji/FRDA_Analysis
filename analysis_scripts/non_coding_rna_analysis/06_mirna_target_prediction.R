suppressPackageStartupMessages({
    library(multiMiR)
    library(dplyr)
    library(tibble)
    library(readr)
})

# 1. LOAD CORRELATION DATA
cor_all <- readRDS("ncrna_analysis/All_ncRNA_mRNA_Correlations.rds")

# 2. FILTER FOR miRNAs AND ANTI-CORRELATION
# Focus on miRNAs that are significantly anti-correlated with potential targets
mir_de <- cor_all |>
    filter(biotype == "miRNA", Correlation < -0.3) |>
    distinct(ncRNA_ID) |>
    pull(ncRNA_ID)

message("Processing targets for ", length(mir_de), " miRNAs...")

# 3. GET PREDICTED TARGETS USING multiMiR
multimir_res <- get_multimir(
    mirna = mir_de,
    summary = TRUE,
    predicted.cutoff = 20, # Top 20% of predicted targets
    predicted.cutoff.type = "p",
    table = "all"
)

targets_df <- multimir_res@data |>
    as_tibble() |>
    select(mature_mirna_id, target_symbol, database, score)

# 4. OVERLAP WITH EMPIRICAL ANTI-CORRELATION
validated_mir_mrna <- cor_all |>
    filter(biotype == "miRNA", Correlation < -0.3) |>
    inner_join(
        targets_df,
        by = c("ncRNA_ID" = "mature_mirna_id", "mRNA_ID" = "target_symbol")
    ) |>
    arrange(Correlation)

# 5. SAVE RESULTS
write.csv(
    validated_mir_mrna,
    "ncrna_analysis/Validated_miRNA_mRNA_Targets.csv",
    row.names = FALSE
)

saveRDS(
    validated_mir_mrna,
    "ncrna_analysis/Validated_miRNA_mRNA_Targets.rds"
)

message("Successfully identified ", nrow(validated_mir_mrna), " high-confidence targets.")