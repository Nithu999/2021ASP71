# ============================================================
# CropStress Intelligence
# Main Workflow Runner
# ============================================================

# Run this script from the project root directory.


# STEP 1: Create required output directories

required_dirs <- c(
  "references",
  "results/tables",
  "results/figures",
  "reports"
)

for (directory in required_dirs) {
  
  if (!dir.exists(directory)) {
    
    dir.create(
      directory,
      recursive = TRUE
    )
    
  }
  
}


# STEP 2: Execute analysis stages in order

cat("\n========================================\n")
cat("Starting CropStress Intelligence Workflow\n")
cat("========================================\n")


# Stage 1: Metadata extraction

cat("\nRunning Stage 1: Metadata Extraction\n")

source("code/01_extract_metadata.R")


# Stage 2: Raw count import

cat("\nRunning Stage 2: Count Data Import\n")

source("code/02_import_counts.R")


# Stage 3: Metadata integration

cat("\nRunning Stage 3: Metadata Integration\n")

source("code/03_integrate_metadata.R")


# Stage 4: Sample quality control

cat("\nRunning Stage 4: Sample Quality Control\n")

source("code/04_sample_quality_control.R")


# Stage 5: Gene identifier inspection

cat("\nRunning Stage 5: Gene Identifier Inspection\n")

source("code/05_gene_identifier_inspection.R")


# Stage 6: Gene annotation investigation

cat("\nRunning Stage 6: Gene Annotation Investigation\n")

source("code/06_gene_annotation_inspection.R")


# Stage 7: Exploratory data analysis

cat("\nRunning Stage 7: Exploratory Data Analysis\n")

source("code/07_exploratory_data_analysis.R")


# Workflow completion

cat("\n========================================\n")
cat("CropStress Intelligence Workflow\n")
cat("All seven stages completed successfully.\n")
cat("========================================\n")