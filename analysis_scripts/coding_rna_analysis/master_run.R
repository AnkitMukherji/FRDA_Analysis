# FRDA Analysis Pipeline - Master Orchestration Script
# Author: Senior Bioinformatician
# Usage: source("analysis_scripts/master_run.R")

message("###################################################")
message("### STARTING FRDA RNA-SEQ ANALYSIS PIPELINE ###")
message("###################################################")

steps <- list(
    "01_data_preprocessing.R",
    "02_batch_correction_pca.R",
    "03_dge_all_batches.R",
    "04_pathway_enrichment.R",
    "05_wgcna.R",
    "06_predictive_modeling.R",
    "07_patient_stratification.R",
    "08_severity_trends.R"
)

for (s in steps) {
    script_path <- file.path("analysis_scripts", s)
    message("\n---------------------------------------------------")
    message(">>> STEP: ", s)
    message("---------------------------------------------------")

    # Track time
    start_time <- Sys.time()
    tryCatch(
        {
            source(script_path)
            end_time <- Sys.time()
            message(">>> SUCCESS: ", s, " (Duration: ", round(difftime(end_time, start_time, units = "mins"), 2), " mins)")
        },
        error = function(e) {
            message(">>> !!! ERROR in ", s, ": ", e$message)
            stop("Pipeline halted due to error in script.")
        }
    )
}

message("\n###################################################")
message("### PIPELINE COMPLETED SUCCESSFULLY ###")
message("###################################################")
