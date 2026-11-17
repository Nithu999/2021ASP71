# ============================================================
# Stage 15: Inspect GFF3 Metadata
# ============================================================

gff_file <- "references/osa1_r7.all_models.gff3.gz"

# Read annotation header

gff_connection <- gzfile(gff_file, open = "rt")

gff_header <- character()

repeat {

  line <- readLines(gff_connection, n = 1)

  if (length(line) == 0) break

  if (!startsWith(line, "#")) break

  gff_header <- c(gff_header, line)

}

close(gff_connection)

cat("\nGFF3 header:\n")
cat(gff_header, sep = "\n")

# Inspect chromosome distribution

cat("\nChromosome distribution:\n")

print(table(gff_table$seqid))

# Inspect gene identifier patterns

cat("\nGene identifier patterns:\n")

print(head(gene_ids, 30))

cat("\nUnique chromosome names:\n")

print(unique(gff_table$seqid))

cat("\nGFF3 metadata inspection completed.\n")