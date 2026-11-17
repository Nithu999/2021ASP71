# ============================================================
# CropStress Intelligence
# Complete Workflow Runner
# ============================================================

# Run from project root.

if (!file.exists("CropStress_Intelligence.Rproj")) {
  stop("Open CropStress_Intelligence.Rproj first.")
}

# SAFETY SWITCH
# Change to TRUE only when you intentionally want to
# regenerate the complete analysis and overwrite outputs.

RUN_FULL_WORKFLOW <- FALSE

if (!RUN_FULL_WORKFLOW) {
  stop(
    "Workflow execution disabled to protect existing results.\n",
    "Review this script and set RUN_FULL_WORKFLOW <- TRUE ",
    "when you intentionally want to rerun everything."
  )
}

# Check required packages
required_packages <- c("readxl", "DESeq2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is missing: ", pkg)
  }
}

# Check required input files
required_files <- c(
  "data/raw/GSE288615_rawCounts_IR64_rehydration_All_samples.txt",
  "references/GSE288615_family.soft.gz",
  "references/12870_2025_7175_MOESM2_ESM.xlsx",
  "references/osa1_r7.all_models.gff3.gz"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# Confirm overwrite risk
cat("\nWARNING: Existing analysis outputs may be overwritten.\n")

confirmation <- readline(
  "Type RUN to execute the complete workflow: "
)

if (confirmation != "RUN") {
  stop("Execution cancelled. Existing results were preserved.")
}

# Required directories
required_dirs <- c(
  "results/tables",
  "results/figures",
  "results/differential_expression",
  "data/processed"
)

for (directory in required_dirs) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# Complete execution order
workflow_stages <- sprintf(
  "code/%02d_%s.R",
  1:21,
  c(
    "extract_metadata",
    "import_counts",
    "integrate_metadata",
    "sample_quality_control",
    "gene_identifier_inspection",
    "gene_annotation_inspection",
    "exploratory_data_analysis",
    "gene_expression_filtering",
    "prepare_differential_expression",
    "validate_differential_expression_inputs",
    "investigate_original_gene_identifiers",
    "investigate_annotation_files",
    "inspect_annotation_compatibility",
    "verify_gene_annotation_order",
    "inspect_gff3_metadata",
    "validate_supplementary_gene_ids",
    "compare_supplementary_and_gff3_ids",
    "data_cleaning_preprocessing",
    "deseq2_differential_expression",
    "validate_differential_expression",
    "outlier_sensitivity_analysis"
  )
)

cat("\nStarting complete CropStress Intelligence workflow.\n")

for (i in seq_along(workflow_stages)) {

  script <- workflow_stages[i]

  if (!file.exists(script)) {
    stop("Missing script: ", script)
  }

  cat("\n========================================\n")
  cat("Stage ", i, " of 21: ", basename(script), "\n")
  cat("========================================\n")

  source(script, local = .GlobalEnv)

  cat("\nCompleted: ", basename(script), "\n")
}

cat("\n========================================\n")
cat("COMPLETE WORKFLOW FINISHED SUCCESSFULLY\n")
cat("========================================\n")

writeLines(
  capture.output(sessionInfo()),
  "reports/workflow_session_info.txt"
)