# ============================================================
# Stage 14: Verify Gene Annotation Order
# ============================================================

# STEP 14.1: Read original count matrix

raw_file <- paste0(
  "data/raw/",
  "GSE288615_rawCounts_IR64_rehydration_All_samples.txt"
)

original_counts <- read.delim(
  raw_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# STEP 14.2: Extract GFF3 gene records

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

colnames(gff_table) <- c(
  "seqid", "source", "type",
  "start", "end", "score",
  "strand", "phase", "attributes"
)

gene_records <- gff_table[
  gff_table$type == "gene",
]

# STEP 14.3: Extract gene IDs

gene_ids <- sub(
  ".*ID=([^;]+).*",
  "\\1",
  gene_records$attributes
)

# STEP 14.4: Validate record counts

cat("\nOriginal count records:\n")
print(nrow(original_counts))

cat("\nGFF3 gene records:\n")
print(length(gene_ids))

# STEP 14.5: Check gene identifier quality

cat("\nAre gene IDs unique?\n")
print(!anyDuplicated(gene_ids))

cat("\nFirst 20 gene IDs:\n")
print(head(gene_ids, 20))

cat("\nLast 20 gene IDs:\n")
print(tail(gene_ids, 20))

# STEP 14.6: Check annotation order

cat("\nFirst 10 annotation coordinates:\n")

print(
  head(
    gene_records[, c("seqid", "start", "end")],
    10
  )
)

cat("\nGene annotation order inspection completed.\n")