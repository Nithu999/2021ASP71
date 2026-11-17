# ============================================================
# Stage 17: Compare Supplementary IDs with GFF3 IDs
# ============================================================

library(readxl)

# STEP 17.1: Read TableS2

excel_file <- "references/12870_2025_7175_MOESM2_ESM.xlsx"

table_s2 <- read_excel(
  excel_file,
  sheet = "TableS2",
  col_names = FALSE
)

# STEP 17.2: Extract supplementary gene IDs

supplementary_ids <- as.character(table_s2[[2]])

supplementary_ids <- supplementary_ids[
  !is.na(supplementary_ids) &
    grepl("^LOC_Os", supplementary_ids)
]

supplementary_ids <- unique(supplementary_ids)

cat("\nSupplementary gene IDs:\n")
print(length(supplementary_ids))

# STEP 17.3: Extract GFF3 gene IDs

gff_file <- "references/osa1_r7.all_models.gff3.gz"

gff_lines <- readLines(gzfile(gff_file))

gff_data <- gff_lines[
  !grepl("^#", gff_lines)
]

gff_table <- read.delim(
  text = gff_data,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

gene_records <- gff_table[
  gff_table$V3 == "gene",
]

gff_ids <- sub(
  ".*ID=([^;]+).*",
  "\\1",
  gene_records$V9
)

# STEP 17.4: Compare identifiers

matching_ids <- intersect(
  supplementary_ids,
  gff_ids
)

cat("\nMatching gene IDs:\n")
print(length(matching_ids))

cat("\nSupplementary IDs absent from GFF3:\n")
print(length(setdiff(supplementary_ids, gff_ids)))

cat("\nGFF3 IDs absent from supplementary table:\n")
print(length(setdiff(gff_ids, supplementary_ids)))

# STEP 17.5: Display examples

cat("\nFirst 10 matching IDs:\n")
print(head(matching_ids, 10))

cat("\nComparison completed.\n")