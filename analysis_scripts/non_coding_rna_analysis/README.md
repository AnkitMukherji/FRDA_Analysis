# Non-coding RNA and regulatory hypotheses

The numbered steps build on prepared expression, metadata, differential-expression tables and co-expression outputs. Review the R environment, expected inputs and path settings before executing any stage. Treat correlations, target predictions and enrichment as evidence for candidate regulatory relationships.

## Script and notebook guide

| File | Purpose |
| --- | --- |
| [`01_data_preparation.R`](01_data_preparation.R) | Prepares expression and metadata inputs for the non-coding RNA stages. |
| [`02_differential_expression.R`](02_differential_expression.R) | Tests differential expression for the configured RNA groups and writes results for follow-up. |
| [`03_linear_modeling.R`](03_linear_modeling.R) | Fits expression models against the selected clinical and technical variables. |
| [`04a_wgcna_analysis.R`](04a_wgcna_analysis.R) | Runs a WGCNA co-expression analysis and exports module information. |
| [`04b_wgcna_analysis.R`](04b_wgcna_analysis.R) | Runs a companion WGCNA analysis with its own analysis settings. Review the two scripts before selecting the appropriate inputs and model. |
| [`05_mrna_lncrna_correlation.R`](05_mrna_lncrna_correlation.R) | Calculates coding/non-coding RNA associations from the prepared expression data. |
| [`06_mirna_target_prediction.R`](06_mirna_target_prediction.R) | Uses anti-correlations to prioritize miRNA pairs and queries multiMiR for target evidence. |
| [`07_functional_enrichment.R`](07_functional_enrichment.R) | Tests functional enrichment on selected gene sets. |
| [`08_network_visualization.R`](08_network_visualization.R) | Builds visual summaries of candidate RNA relationships. |
| [`09_wgcna_functional_enrichment.R`](09_wgcna_functional_enrichment.R) | Evaluates functional enrichment for genes in the WGCNA modules. |

[Repository overview](../../README.md)
