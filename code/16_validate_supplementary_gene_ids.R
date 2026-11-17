# ============================================================
# Stage 16: Validate Supplementary Gene Identifiers
# ============================================================

library(readxl)

# STEP 16.1: Read TableS2 without assuming column headers

excel_file <- "references/12870_2025_7175_MOESM2_ESM.xlsx"

table_s2 <- read_excel(
  excel_file,
  sheet = "TableS2",
  col_names = FALSE
)

cat("\nTableS2 dimensions:\n")
print(dim(table_s2))

cat("\nFirst 20 rows of TableS2:\n")
print(table_s2[1:min(20, nrow(table_s2)), ])

# STEP 16.2: Extract gene IDs from GFF3

gff_file <- "references/osa1_r7.all_models.gff3.gz"

gff_lines <- readLines(gzfile(gff_file))

gff_data <- gff_lines[!grepl("^#", gff_lines)]

gff_table <- read.delim(
  text = gff_data,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

gene_records <- gff_table[gff_table$V3 == "gene", ]

gff_gene_ids <- sub(
  ".*ID=([^;]+).*",
  "\\1",
  gene_records$V9
)

cat("\nTotal GFF3 gene IDs:\n")
print(length(gff_gene_ids))

cat("\nSupplementary identifier inspection completed.\n")