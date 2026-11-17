# ============================================================
# Stage 13: Inspect Annotation Compatibility
# ============================================================

# STEP 13.1: Inspect GFF3 annotation records

gff_file <- "references/osa1_r7.all_models.gff3.gz"

gff <- readLines(gzfile(gff_file))

# Remove comment lines

gff_data <- gff[!grepl("^#", gff)]

# Split GFF3 columns

gff_table <- read.delim(
  text = gff_data,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(gff_table) <- c(
  "seqid", "source", "type",
  "start", "end", "score",
  "strand", "phase", "attributes"
)

# STEP 13.2: Count annotation feature types

cat("\nGFF3 feature types:\n")
print(table(gff_table$type))

# STEP 13.3: Extract gene records

gene_records <- gff_table[
  gff_table$type == "gene",
]

cat("\nNumber of gene records:\n")
print(nrow(gene_records))

# STEP 13.4: Inspect gene identifiers

cat("\nFirst 10 gene annotation records:\n")

print(head(
  gene_records[, c("seqid", "start", "end", "attributes")],
  10
))

# STEP 13.5: Compare record counts

cat("\nCount matrix records:\n")
print(nrow(filtered_counts))

cat("\nOriginal count matrix records:\n")
print(nrow(counts))

cat("\nGFF3 gene records:\n")
print(nrow(gene_records))

cat("\nAnnotation compatibility inspection completed.\n")