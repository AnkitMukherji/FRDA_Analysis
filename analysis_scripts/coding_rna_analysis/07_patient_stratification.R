# FRDA Analysis Pipeline - 07 Patient Stratification
# Author: Senior Bioinformatician

# --- Load Config & Utils ---
source("analysis_scripts/00_init.R")
source("analysis_scripts/utils.R")

# --- Libraries ---
library(tidyverse)
library(ggpubr)
library(rstatix)

# --- Load Data ---
data_corrected <- readRDS(file.path(DIR_RES, "data_corrected.rds"))
meta <- data_corrected$meta

# --- 1. Composite Severity Index ---
message(">>> Calculating Composite Severity Index...")
# Focus on Patients with Homozygous GAA1
patients_hom <- meta |>
    filter(Condition == "Patient") |>
    # Heuristic for homozygous based on control repeat length (from .qmd line 2297)
    filter(as.numeric(GAA1) > 40)

severity_df <- patients_hom |>
    select(`FSA scale`, `mFARS total`, `mFARS USS`) |>
    mutate(across(everything(), as.numeric))

# Z-score each scale
z_mat <- as.data.frame(scale(severity_df))
colnames(z_mat) <- paste0(colnames(z_mat), "_z")

severity_df <- bind_cols(severity_df, z_mat) |>
    mutate(
        n_present = rowSums(!is.na(across(ends_with("_z")))),
        index = rowMeans(across(ends_with("_z")), na.rm = TRUE)
    )

# --- 2. Tertile Grouping ---
message(">>> Assigning Tertile Severity Groups...")
# Use samples with at least 2 scales
eligible <- severity_df |> filter(n_present >= 2)
cutoffs <- quantile(eligible$index, probs = c(1 / 3, 2 / 3), na.rm = TRUE)

severity_df <- severity_df |>
    mutate(severity_group = case_when(
        n_present < 2 ~ "Insufficient data",
        index <= cutoffs[1] ~ "Mild",
        index <= cutoffs[2] ~ "Moderate",
        TRUE ~ "Severe"
    )) |>
    mutate(severity_group = factor(severity_group, levels = c("Mild", "Moderate", "Severe", "Insufficient data")))

# Merge back to meta
meta_strat <- meta
meta_strat$severity_group <- severity_df$severity_group[match(rownames(meta_strat), rownames(severity_df))]
meta_strat$severity_index <- severity_df$index[match(rownames(meta_strat), rownames(severity_df))]

# --- 3. Visualization ---
message(">>> Visualizing Stratification Outcomes...")
p_index <- ggplot(
    severity_df |> filter(severity_group != "Insufficient data"),
    aes(x = index, fill = severity_group)
) +
    geom_histogram(bins = 20, alpha = 0.7, color = "black") +
    scale_fill_manual(values = PAL_SEVERITY) +
    theme_bw() +
    labs(title = "Composite Severity Index Distribution", x = "Mean Z-score Index")

ggsave(file.path(DIR_ML, "severity_index_distribution.png"), p_index, width = 10, height = 7)

saveRDS(meta_strat, file.path(DIR_RES, "meta_stratified.rds"))

message(">>> Patient Stratification Complete.")
