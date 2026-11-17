# ============================================================
# CropStress Intelligence
# Stage 6: Gene Annotation Investigation
# ============================================================

# STEP 1: Define input and output files

annotation_file <- file.path(
  "references",
  "12870_2025_7175_MOESM2_ESM.xlsx"
)

geo_metadata_file <- file.path(
  "references",
  "GSE288615_family.soft.gz"
)

gff_file <- file.path(
  "references",
  "osa1_r7.all_models.gff3.gz"
)

counts_file <- file.path(
  "data",
  "raw",
  "GSE288615_rawCounts_IR64_rehydration_All_samples.txt"
)

output_directory <- file.path(
  "results",
  "tables"
)

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Check required input files

required_files <- c(
  annotation_file,
  geo_metadata_file,
  gff_file,
  counts_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Required files are missing:",
      paste(missing_files, collapse = "\n")
    )
  )
}

# Load required package

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install readxl using install.packages('readxl').")
}

library(readxl)


# ============================================================
# STEP 2: Inspect supplementary annotation workbook
# ============================================================

sheet_names <- excel_sheets(annotation_file)

sheet_information <- data.frame(
  sheet_name = character(),
  rows = integer(),
  columns = integer()
)

for (sheet in sheet_names) {
  
  sheet_data <- read_excel(
    annotation_file,
    sheet = sheet
  )
  
  sheet_information <- rbind(
    sheet_information,
    data.frame(
      sheet_name = sheet,
      rows = nrow(sheet_data),
      columns = ncol(sheet_data)
    )
  )
}

write.csv(
  sheet_information,
  file.path(
    output_directory,
    "gene_annotation_workbook_summary.csv"
  ),
  row.names = FALSE
)

cat("\nSupplementary workbook information:\n")
print(sheet_information)


# ============================================================
# STEP 3: Extract supplementary gene identifiers from TableS2
# ============================================================

if (!"TableS2" %in% sheet_names) {
  stop("TableS2 was not found in the supplementary workbook.")
}

table_s2 <- read_excel(
  annotation_file,
  sheet = "TableS2"
)

cat("\nTableS2 dimensions:\n")
print(dim(table_s2))

cat("\nTableS2 column names:\n")
print(names(table_s2))

# The previous investigation identified gene IDs
# in the second column of TableS2.

if (ncol(table_s2) < 2) {
  stop("TableS2 does not contain the expected second column.")
}

gene_ids_s2_all <- as.character(table_s2[[2]])

# Extract valid rice gene identifiers

supplementary_gene_ids <- unique(
  gene_ids_s2_all[
    !is.na(gene_ids_s2_all) &
      grepl(
        "^LOC_Os[0-9]+g[0-9]+$",
        gene_ids_s2_all
      )
  ]
)

cat("\nSupplementary gene identifiers:\n")
print(length(supplementary_gene_ids))

cat("\nFirst 15 supplementary gene IDs:\n")
print(head(supplementary_gene_ids, 15))

cat("\nDuplicated supplementary gene IDs:\n")
print(sum(duplicated(supplementary_gene_ids)))

write.csv(
  data.frame(
    gene_id = supplementary_gene_ids
  ),
  file.path(
    output_directory,
    "supplementary_gene_identifiers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# STEP 4: Investigate GEO processing information
# ============================================================

geo_connection <- gzfile(
  geo_metadata_file,
  open = "rt"
)

soft_lines <- readLines(
  geo_connection,
  warn = FALSE
)

close(geo_connection)

processing_lines <- grep(
  "data_processing|genome|assembly|featureCounts|STAR|supplementary",
  soft_lines,
  value = TRUE,
  ignore.case = TRUE
)

cat("\nGEO processing information:\n")

cat(
  processing_lines,
  sep = "\n"
)

writeLines(
  processing_lines,
  file.path(
    output_directory,
    "geo_processing_information.txt"
  )
)


# ============================================================
# STEP 5: Read and inspect MSU7 GFF3 annotation
# ============================================================

gff_connection <- gzfile(
  gff_file,
  open = "rt"
)

gff_lines <- readLines(
  gff_connection,
  warn = FALSE
)

close(gff_connection)

# Remove comment and empty lines

gff_data <- gff_lines[
  !grepl("^#", gff_lines) &
    trimws(gff_lines) != ""
]

# Split GFF3 records into columns

gff_columns <- strsplit(
  gff_data,
  "\t",
  fixed = TRUE
)

valid_records <- lengths(gff_columns) == 9

gff_data <- gff_data[valid_records]
gff_columns <- gff_columns[valid_records]

# Extract feature types

feature_types <- vapply(
  gff_columns,
  function(x) x[3],
  character(1)
)

cat("\nGFF3 feature distribution:\n")
print(sort(table(feature_types), decreasing = TRUE))

cat("\nTotal annotation records:\n")
print(length(gff_data))


# ============================================================
# STEP 6: Extract MSU7 gene identifiers
# ============================================================

gene_records <- gff_data[
  feature_types == "gene"
]

gene_columns <- gff_columns[
  feature_types == "gene"
]

# Extract ID attribute from GFF3

gene_attributes <- vapply(
  gene_columns,
  function(x) x[9],
  character(1)
)

extract_gff_id <- function(attributes) {
  
  id_match <- regmatches(
    attributes,
    regexpr(
      "(^|;)ID=[^;]+",
      attributes
    )
  )
  
  if (length(id_match) == 0 || id_match == "") {
    return(NA_character_)
  }
  
  sub(
    "(^|;)ID=",
    "",
    id_match
  )
}

gene_ids_msu7 <- vapply(
  gene_attributes,
  extract_gff_id,
  character(1)
)

# Validate extracted identifiers

cat("\nMSU7 gene annotation summary:\n")

cat("\nTotal gene records:\n")
print(length(gene_records))

cat("\nExtracted gene identifiers:\n")
print(length(gene_ids_msu7))

cat("\nMissing gene identifiers:\n")
print(sum(is.na(gene_ids_msu7)))

cat("\nEmpty gene identifiers:\n")
print(sum(gene_ids_msu7 == "", na.rm = TRUE))

cat("\nDuplicated gene identifiers:\n")
print(sum(duplicated(gene_ids_msu7)))

cat("\nFirst 15 gene identifiers:\n")
print(head(gene_ids_msu7, 15))

cat("\nLast 15 gene identifiers:\n")
print(tail(gene_ids_msu7, 15))


# ============================================================
# STEP 7: Investigate gene identifier categories
# ============================================================

loc_os_ids <- gene_ids_msu7[
  grepl(
    "^LOC_Os",
    gene_ids_msu7
  )
]

other_gene_ids <- gene_ids_msu7[
  !grepl(
    "^LOC_Os",
    gene_ids_msu7
  ) &
    !is.na(gene_ids_msu7)
]

cat("\nLOC_Os gene identifiers:\n")
print(length(loc_os_ids))

cat("\nOther gene identifiers:\n")
print(length(other_gene_ids))

cat("\nFirst 20 other identifiers:\n")
print(head(other_gene_ids, 20))

# Save complete MSU7 gene annotation identifiers

msu7_annotation <- data.frame(
  annotation_row = seq_along(gene_ids_msu7),
  gene_id = gene_ids_msu7,
  stringsAsFactors = FALSE
)

write.csv(
  msu7_annotation,
  file.path(
    output_directory,
    "msu7_gene_identifiers.csv"
  ),
  row.names = FALSE
)


# ============================================================
# STEP 8: Validate gene identifier extraction
# ============================================================

valid_gene_ids <- !is.na(gene_ids_msu7) &
  gene_ids_msu7 != ""

identifier_summary <- data.frame(
  total_gene_records = length(gene_ids_msu7),
  valid_identifiers = sum(valid_gene_ids),
  missing_identifiers = sum(is.na(gene_ids_msu7)),
  empty_identifiers = sum(gene_ids_msu7 == "", na.rm = TRUE),
  unique_identifiers = length(unique(gene_ids_msu7[valid_gene_ids])),
  duplicate_identifiers = sum(duplicated(gene_ids_msu7[valid_gene_ids]))
)

cat("\nGene identifier validation:\n")
print(identifier_summary)

write.csv(
  identifier_summary,
  file.path(
    output_directory,
    "gene_identifier_validation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# STEP 9: Compare supplementary IDs with MSU7
# ============================================================

matching_gene_ids <- intersect(
  supplementary_gene_ids,
  gene_ids_msu7[valid_gene_ids]
)

unmatched_gene_ids <- setdiff(
  supplementary_gene_ids,
  gene_ids_msu7[valid_gene_ids]
)

cat("\nSupplementary unique gene IDs:\n")
print(length(supplementary_gene_ids))

cat("\nIDs matching MSU7:\n")
print(length(matching_gene_ids))

cat("\nIDs not matching MSU7:\n")
print(length(unmatched_gene_ids))

cat("\nFirst 20 unmatched IDs:\n")
print(head(unmatched_gene_ids, 20))

# Calculate matching percentage

matching_percentage <- (
  length(matching_gene_ids) /
    length(supplementary_gene_ids)
) * 100

cat("\nSupplementary identifier matching percentage:\n")
print(round(matching_percentage, 2))

# Save comparison results

identifier_comparison <- data.frame(
  total_supplementary_ids = length(supplementary_gene_ids),
  matching_ids = length(matching_gene_ids),
  unmatched_ids = length(unmatched_gene_ids),
  matching_percentage = round(matching_percentage, 2)
)

write.csv(
  identifier_comparison,
  file.path(
    output_directory,
    "supplementary_msu7_identifier_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    gene_id = unmatched_gene_ids
  ),
  file.path(
    output_directory,
    "unmatched_supplementary_gene_ids.csv"
  ),
  row.names = FALSE
)


# ============================================================
# STEP 10: Inspect original RNA-seq count matrix
# ============================================================

counts <- read.delim(
  counts_file,
  header = TRUE,
  check.names = FALSE
)

cat("\nOriginal count matrix dimensions:\n")
print(dim(counts))

cat("\nNumber of samples:\n")
print(ncol(counts))

cat("\nNumber of expression records:\n")
print(nrow(counts))

cat("\nFirst five sample names:\n")
print(head(colnames(counts), 5))


# ============================================================
# STEP 11: Compare count matrix and annotation dimensions
# ============================================================

count_matrix_genes <- nrow(counts)

annotation_genes <- length(gene_ids_msu7)

cat("\nCount matrix records:\n")
print(count_matrix_genes)

cat("\nMSU7 gene annotation records:\n")
print(annotation_genes)

dimension_match <- count_matrix_genes == annotation_genes

cat("\nDo the dimensions match?\n")
print(dimension_match)

# Create a candidate mapping only for investigation.
# Matching dimensions alone do not establish matching gene order.

if (dimension_match &&
    all(valid_gene_ids) &&
    !anyDuplicated(gene_ids_msu7)) {
  
  candidate_mapping <- data.frame(
    count_row = seq_len(count_matrix_genes),
    gene_id = gene_ids_msu7,
    mapping_status = "Order-based candidate; requires verification",
    stringsAsFactors = FALSE
  )
  
  write.csv(
    candidate_mapping,
    file.path(
      output_directory,
      "candidate_count_gene_mapping.csv"
    ),
    row.names = FALSE
  )
  
  cat(
    "\nCandidate mapping saved for further verification.\n"
  )
  
} else {
  
  cat(
    "\nAutomatic candidate mapping was not created.\n"
  )
  
}


# ============================================================
# STEP 12: Save final annotation investigation summary
# ============================================================

final_summary <- data.frame(
  dataset = "GSE288615",
  organism = "Oryza sativa",
  cultivar = "IR64",
  count_matrix_genes = count_matrix_genes,
  count_matrix_samples = ncol(counts),
  msu7_gene_records = annotation_genes,
  valid_msu7_gene_ids = sum(valid_gene_ids),
  unique_msu7_gene_ids = length(unique(gene_ids_msu7[valid_gene_ids])),
  supplementary_gene_ids = length(supplementary_gene_ids),
  matching_supplementary_ids = length(matching_gene_ids),
  unmatched_supplementary_ids = length(unmatched_gene_ids),
  matching_percentage = round(matching_percentage, 2),
  count_annotation_dimensions_match = dimension_match
)

write.csv(
  final_summary,
  file.path(
    output_directory,
    "gene_annotation_investigation_summary.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("Gene annotation investigation completed.\n")
cat("============================================\n")

cat("\nFinal summary:\n")
print(final_summary)

cat("\nOutput files saved in results/tables/.\n")