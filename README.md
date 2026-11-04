# Friedreich Ataxia: Blood Transcriptomics and Regulatory Networks

Research analysis of the blood transcriptomic landscape of **Friedreich ataxia (FRDA)**, connecting expression changes with pathway activity, co-expression networks, and clinical severity.

The repository contains R scripts for coding and non-coding RNA analyses. It complements the [cardiomyocyte project](https://github.com/AnkitMukherji/FRDA_Cardiomyopathy_Analysis), which examines cardiac involvement and differentiation stage.

## Scientific questions

- Which transcriptional changes distinguish patients and controls after accounting for technical variation?
- Which pathways and co-expression modules track disease status or severity?
- How do mRNA, lncRNA, and miRNA analyses generate complementary regulatory hypotheses?

## Analysis map

| Location | Analysis |
| --- | --- |
| [Coding RNA scripts](analysis_scripts/coding_rna_analysis/) | Preprocessing, batch assessment/PCA, differential expression, pathway enrichment, WGCNA, exploratory predictive modeling, stratification, and severity trends |
| [Non-coding RNA scripts](analysis_scripts/non_coding_rna_analysis/) | Data preparation, differential expression, linear models, WGCNA, mRNA–lncRNA correlations, miRNA target prediction, enrichment, and network visualization |
| [Shared coding RNA helpers](analysis_scripts/coding_rna_analysis/00_helper_scripts.R) | Metadata handling, model-design checks, plotting, differential expression, and pathway utilities |
| [renv.lock](renv.lock) / [Dockerfile](Dockerfile) | Recorded R dependencies and container definition |

## Inputs and reproducibility

These are study-specific analysis scripts, with external expression objects, metadata workbooks, clinical severity tables, and QC summaries. The preprocessing script reads Salmon-related RDS objects from `salmon/`, workbooks from `metadata/`, and QC information from `QC/`.

Restore the recorded R packages from the repository root:

```r
renv::restore()
```

**The current checkout requires configuration before execution.** The preprocessing script sources `helper_scripts.R`, while the committed helper is named `00_helper_scripts.R`. Later scripts also reference `analysis_scripts/00_init.R` and `analysis_scripts/utils.R`, which are absent from this checkout. Resolve these references, supply study inputs, and review output directories before running stages individually. The container definition does not establish that the full analysis runs on a fresh machine.

## Interpreting the analyses

Assess effect sizes, uncertainty, batch structure, and covariate balance together. Correlated expression and pathway enrichment generate biological hypotheses; they do not establish regulatory direction or causality.

The predictive-modeling script preselects signature genes using differential-expression results before its cross-validation loop. Performance from that workflow is exploratory: unbiased evaluation requires feature selection within each training fold and independent validation. No diagnostic performance claim is made here.

## Research context

Ankit Mukherjee · computational biology and disease transcriptomics. See the [research portfolio](https://github.com/AnkitMukherji/CV) for experience and related work.
