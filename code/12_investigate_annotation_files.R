# ============================================================
# Stage 12: Investigate Annotation Files
# ============================================================

# STEP 12.1: Inspect supplementary Excel workbook

library(readxl)

excel_file <- "references/12870_2025_7175_MOESM2_ESM.xlsx"

cat("\nExcel workbook sheets:\n")

sheets <- excel_sheets(excel_file)

print(sheets)

# STEP 12.2: Inspect dimensions and column names

for (sheet in sheets) {

  cat("\n====================================\n")
  cat("Sheet:", sheet, "\n")

  data <- read_excel(
    excel_file,
    sheet = sheet,
    n_max = 5
  )

  cat("\nDimensions of preview:\n")
  print(dim(data))

  cat("\nColumn names:\n")
  print(names(data))

  cat("\nFirst five rows:\n")
  print(data)

}

# STEP 12.3: Inspect GFF3 annotation file

gff_file <- "references/osa1_r7.all_models.gff3.gz"

cat("\nGFF3 file exists:\n")
print(file.exists(gff_file))

gff_preview <- readLines(
  gzfile(gff_file),
  n = 15
)

cat("\nFirst 15 GFF3 lines:\n")

cat(gff_preview, sep = "\n")

cat("\nAnnotation investigation completed.\n")