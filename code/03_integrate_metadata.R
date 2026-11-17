# CropStress Intelligence
# Stage 3: Map GEO sample accessions to count-matrix library names
# Output: results/tables/integrated_sample_metadata.csv

# STEP 14: Build official GEO-to-count matrix mapping


library_names <- vapply(
  seq_along(sample_starts),
  function(i) {

    start <- sample_starts[i]

    end <- if (i < length(sample_starts)) {
      sample_starts[i + 1] - 1
    } else {
      length(soft_lines)
    }

    sample_block <- soft_lines[start:end]

    library_line <- grep(
      "^!Sample_description = Library name:",
      sample_block,
      value = TRUE
    )

    if (length(library_line) == 0) {
      return(NA_character_)
    }

    sub(
      "^!Sample_description = Library name:\\s*",
      "",
      library_line[1]
    )
  },
  character(1)
)


sample_mapping <- data.frame(

  geo_accession = sample_metadata$geo_accession,

  library_name = library_names,

  stringsAsFactors = FALSE
)


cat("\nFirst ten GEO-to-library mappings:\n")

print(head(sample_mapping, 10))


# STEP 15: Validate sample mapping


cat(
  "\nTotal mapping records:",
  nrow(sample_mapping),
  "\n"
)


cat(
  "Missing library names:",
  sum(is.na(sample_mapping$library_name)),
  "\n"
)


cat(
  "Duplicate library names:",
  sum(duplicated(sample_mapping$library_name)),
  "\n"
)


cat(
  "All libraries found in count matrix:",
  all(sample_mapping$library_name %in% colnames(counts)),
  "\n"
)


cat(
  "Complete sample identifier match:",
  setequal(
    sample_mapping$library_name,
    colnames(counts)
  ),
  "\n"
)


cat(
  "Sample order identical:",
  identical(
    sample_mapping$library_name,
    colnames(counts)
  ),
  "\n"
)


cat("\nLibraries missing from count matrix:\n")

print(
  setdiff(
    sample_mapping$library_name,
    colnames(counts)
  )
)

cat("\nCount matrix samples missing from GEO mapping:\n")

print(
  setdiff(
    colnames(counts),
    sample_mapping$library_name
  )
)


# STEP 16: Create integrated sample annotation table


integrated_metadata <- sample_metadata

integrated_metadata$library_name <-
  sample_mapping$library_name


stopifnot(
  identical(
    integrated_metadata$library_name,
    colnames(counts)
  )
)


integrated_metadata <- integrated_metadata[
  ,
  c(
    "geo_accession",
    "library_name",
    "cultivar",
    "tissue",
    "treatment",
    "day",
    "replicate",
    "title"
  )
]


cat("\nIntegrated sample annotation table:\n")

print(head(integrated_metadata, 10))


cat("\nIntegrated metadata dimensions:\n")

print(dim(integrated_metadata))


cat(
  "\nSample mapping verified:",
  identical(
    integrated_metadata$library_name,
    colnames(counts)
  ),
  "\n"
)


# STEP 17: Save integrated sample annotation table


output_directory <- file.path(
  "results",
  "tables"
)


if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}


metadata_output <- file.path(
  output_directory,
  "integrated_sample_metadata.csv"
)


write.csv(
  integrated_metadata,
  file = metadata_output,
  row.names = FALSE
)


cat("\nSaved metadata file:\n")
print(metadata_output)

cat("\nFile size in bytes:\n")
print(file.info(metadata_output)$size)


saved_metadata <- read.csv(
  metadata_output,
  stringsAsFactors = FALSE
)

cat("\nSaved metadata dimensions:\n")
print(dim(saved_metadata))

