# Coding RNA in Friedreich ataxia

Start with the preprocessing script and follow the numbered stages as their inputs become available. Use the repository's R environment and inspect the root README for missing helper paths and external study inputs. Most stages consume intermediate RDS objects or analysis tables.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`00_helper_scripts.R`](00_helper_scripts.R) | Shared utilities for metadata, design checks, expression models, enrichment and visualization. |
| [`01_data_preprocessing.R`](01_data_preprocessing.R) | Combines Salmon-related expression objects with sample and severity metadata and applies expression and sample QC filters before saving RDS inputs. |
| [`02_batch_correction_pca.R`](02_batch_correction_pca.R) | Transforms counts, evaluates sequencing batch and disease structure with PCA, and saves the adjusted analysis object. |
| [`03_dge_all_batches.R`](03_dge_all_batches.R) | Runs differential expression across prepared batches and saves gene-level results for follow-up. |
| [`04_pathway_enrichment.R`](04_pathway_enrichment.R) | Uses differential expression lists and rankings for pathway enrichment and related plots. |
| [`05_wgcna.R`](05_wgcna.R) | Builds weighted gene co-expression modules and relates module summaries to sample traits. |
| [`06_predictive_modeling.R`](06_predictive_modeling.R) | Fits exploratory classifiers to expression features and examines classification performance. Check how feature selection interacts with validation splits before interpreting accuracy. |
| [`07_patient_stratification.R`](07_patient_stratification.R) | Summarizes expression-based patient stratification within the specified GAA repeat subgroup and compares its clinical characteristics. |
| [`08_severity_trends.R`](08_severity_trends.R) | Models expression trends across clinical severity using the configured covariates and contrasts. |
| [`paper_plots.R`](paper_plots.R) | Assembles publication-oriented figures from expression matrices and saved model results. |

[Repository overview](../../README.md)
