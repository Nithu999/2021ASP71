
# ============================================================
# CropStress Intelligence 🌾
# Script: 01_dataset_metadata.R
# Purpose: Download, extract and validate GEO sample metadata
# Dataset: GSE288615
# ============================================================


# ============================================================
# STEP 1: Define official GEO download URL
# ============================================================

geo_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE288nnn/GSE288615/soft/",
  "GSE288615_family.soft.gz"
)


# ============================================================
# STEP 2: Define metadata output location
# ============================================================
metadata_file <- file.path(
  "references",
  "GSE288615_family.soft.gz"
)


# ============================================================
# STEP 3: Check references directory
# ============================================================

if (!dir.exists("references")) {
  dir.create("references")
}


# ============================================================
# STEP 4: Download official GEO metadata
# ============================================================

download.file(
  url = geo_url,
  destfile = metadata_file,
  mode = "wb"
)


# ============================================================
# STEP 5: Verify downloaded file
# ============================================================

cat("\nDownloaded metadata file size:\n")

print(file.info(metadata_file)$size)


# ============================================================
# STEP 6: Read compressed GEO metadata
# ============================================================

geo_connection <- gzfile(
  metadata_file,
  open = "rt"
)

soft_lines <- readLines(
  geo_connection,
  warn = FALSE
)

close(geo_connection)


# ============================================================
# STEP 7: Count official GEO sample records
# ============================================================

sample_records <- grep(
  "^\\^SAMPLE =",
  soft_lines,
  value = TRUE
)

cat(
  "\nNumber of GEO sample records:",
  length(sample_records),
  "\n"
)


# ============================================================
# STEP 8: Inspect sample annotation fields
# ============================================================

sample_annotations <- grep(
  "^!Sample_(geo_accession|title|source_name_ch1|characteristics_ch1)",
  soft_lines,
  value = TRUE
)

cat("\nFirst 30 sample annotation entries:\n")

cat(
  head(sample_annotations, 30),
  sep = "\n"
)


# ============================================================
# STEP 9: Extract individual GEO sample records
# ============================================================

# Find the beginning of every sample record

sample_starts <- grep(
  "^\\^SAMPLE =",
  soft_lines
)


# Function to extract a specific sample field

extract_field <- function(lines, pattern) {
  
  result <- grep(
    pattern,
    lines,
    value = TRUE
  )
  
  if (length(result) == 0) {
    return(NA_character_)
  }
  
  sub(
    "^[^=]*=\\s*",
    "",
    result[1]
  )
}


# Process all GEO sample records

sample_metadata_list <- lapply(
  seq_along(sample_starts),
  function(i) {
    
    start <- sample_starts[i]
    
    end <- if (i < length(sample_starts)) {
      sample_starts[i + 1] - 1
    } else {
      length(soft_lines)
    }
    
    sample_block <- soft_lines[start:end]
    
    
    # Extract sample characteristics
    
    characteristics <- grep(
      "^!Sample_characteristics_ch1",
      sample_block,
      value = TRUE
    )
    
    
    # Function to extract individual characteristics
    
    get_characteristic <- function(key) {
      
      pattern <- paste0(
        "^!Sample_characteristics_ch1 = ",
        key,
        ":"
      )
      
      result <- grep(
        pattern,
        characteristics,
        value = TRUE
      )
      
      if (length(result) == 0) {
        return(NA_character_)
      }
      
      trimws(
        sub(
          paste0("^.*", key, ":\\s*"),
          "",
          result[1]
        )
      )
    }
    
    
    # Create one structured record per sample
    
    data.frame(
      
      geo_accession = extract_field(
        sample_block,
        "^!Sample_geo_accession"
      ),
      
      title = extract_field(
        sample_block,
        "^!Sample_title"
      ),
      
      tissue = get_characteristic("tissue"),
      
      cultivar = get_characteristic("cultivar"),
      
      treatment = get_characteristic("treatment"),
      
      day = get_characteristic("days after_treatment"),
      
      replicate = get_characteristic("replicate"),
      
      stringsAsFactors = FALSE
    )
    
  }
)


# Combine all sample records into one data frame

sample_metadata <- do.call(
  rbind,
  sample_metadata_list
)


# Display first six sample records

cat("\nFirst six sample records:\n")

print(head(sample_metadata))


# Check metadata dimensions

cat("\nMetadata dimensions:\n")

print(dim(sample_metadata))


# ============================================================
# STEP 10: Validate extracted sample metadata
# ============================================================

cat("\nMetadata validation:\n")


# Check total sample records

cat(
  "Total samples:",
  nrow(sample_metadata),
  "\n"
)


# Check missing values

cat(
  "Missing values:",
  sum(is.na(sample_metadata)),
  "\n"
)


# Check duplicate GEO accessions

cat(
  "Duplicate GEO accessions:",
  sum(duplicated(sample_metadata$geo_accession)),
  "\n"
)


# Check unique GEO accessions

cat(
  "Unique GEO accessions:",
  length(unique(sample_metadata$geo_accession)),
  "\n"
)


# ============================================================
# STEP 11: Summarise experimental design
# ============================================================


# 1. Identify treatment groups

cat("\nTreatment groups:\n")

print(unique(sample_metadata$treatment))


# 2. Identify collection time points
# Preserve original labels such as "0 (2h)"

cat("\nCollection time points:\n")

print(unique(sample_metadata$day))


# 3. Identify tissue types

cat("\nTissue types:\n")

print(unique(sample_metadata$tissue))


# 4. Count samples per treatment

cat("\nSamples per treatment:\n")

print(table(sample_metadata$treatment))


# 5. Examine treatment, day and tissue combinations

cat("\nExperimental group distribution:\n")

group_distribution <- xtabs(
  ~ treatment + day + tissue,
  data = sample_metadata
)

print(group_distribution)



# ============================================================
# STEP 12: Import and inspect raw RNA-seq count matrix
# ============================================================

# 1. Define raw count file location

counts_file <- file.path(
  "data",
  "raw",
  "GSE288615_rawCounts_IR64_rehydration_All_samples.txt"
)

# 2. Import count matrix
# Do not assign row.names yet.
# We first need to understand the file structure.

counts <- read.delim(
  counts_file,
  header = TRUE,
  check.names = FALSE
)

# 3. Check matrix dimensions

cat("\nCount matrix dimensions:\n")
print(dim(counts))

# 4. Inspect column names

cat("\nFirst 10 column names:\n")
print(head(colnames(counts), 10))

# 5. Inspect first six rows

cat("\nFirst six rows:\n")
print(head(counts))

# 6. Inspect data structure

cat("\nCount matrix structure:\n")
str(counts)

# 7. Check whether all columns contain numeric values

cat("\nAre all columns numeric?\n")

print(
  all(sapply(counts, is.numeric))
)


# 8. Verify expression matrix dimensions

cat("\nExpression matrix dimensions:\n")
print(dim(counts))

# 9. Count samples

cat("\nNumber of samples:\n")
print(ncol(counts))

# 10. Count expression records

cat("\nNumber of expression records:\n")
print(nrow(counts))

# 11. Display all sample identifiers

cat("\nSample identifiers:\n")
print(colnames(counts))

# 12. Inspect first five expression records

cat("\nFirst five expression records:\n")
print(head(counts, 5))


# ============================================================
# STEP 13: Search GEO metadata for sample mapping information
# ============================================================

# Search for library names and count matrix column names

mapping_lines <- grep(
  "Library name|Column name in rawCounts|IR0C1L1",
  soft_lines,
  value = TRUE,
  ignore.case = TRUE
)

# Display matching metadata entries

cat("\nOfficial sample mapping information:\n")

print(head(mapping_lines, 30))

# Count matching entries

cat(
  "\nTotal matching entries:",
  length(mapping_lines),
  "\n"
)


# ============================================================
# STEP 14: Build official GEO-to-count matrix mapping
# ============================================================

# Extract library names from individual GEO sample records

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


# Create the official sample mapping table

sample_mapping <- data.frame(
  
  geo_accession = sample_metadata$geo_accession,
  
  library_name = library_names,
  
  stringsAsFactors = FALSE
)


# Display first ten mapping records

cat("\nFirst ten GEO-to-library mappings:\n")

print(head(sample_mapping, 10))


# ============================================================
# STEP 15: Validate sample mapping
# ============================================================

# Check total mapping records

cat(
  "\nTotal mapping records:",
  nrow(sample_mapping),
  "\n"
)


# Check missing library names

cat(
  "Missing library names:",
  sum(is.na(sample_mapping$library_name)),
  "\n"
)


# Check duplicate library names

cat(
  "Duplicate library names:",
  sum(duplicated(sample_mapping$library_name)),
  "\n"
)


# Check whether all library names exist in count matrix

cat(
  "All libraries found in count matrix:",
  all(sample_mapping$library_name %in% colnames(counts)),
  "\n"
)


# Check whether both datasets contain exactly the same identifiers

cat(
  "Complete sample identifier match:",
  setequal(
    sample_mapping$library_name,
    colnames(counts)
  ),
  "\n"
)


# Check whether the sample order is identical

cat(
  "Sample order identical:",
  identical(
    sample_mapping$library_name,
    colnames(counts)
  ),
  "\n"
)


# Identify any unmatched identifiers

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


# ============================================================
# STEP 16: Create integrated sample annotation table
# ============================================================

# Add official library names to GEO sample metadata

integrated_metadata <- sample_metadata

integrated_metadata$library_name <-
  sample_mapping$library_name

# Verify library names match count matrix columns

stopifnot(
  identical(
    integrated_metadata$library_name,
    colnames(counts)
  )
)

# Arrange columns in a meaningful order

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

# Display integrated metadata

cat("\nIntegrated sample annotation table:\n")

print(head(integrated_metadata, 10))

# Verify dimensions

cat("\nIntegrated metadata dimensions:\n")

print(dim(integrated_metadata))

# Verify sample matching

cat(
  "\nSample mapping verified:",
  identical(
    integrated_metadata$library_name,
    colnames(counts)
  ),
  "\n"
)


# ============================================================
# STEP 17: Save integrated sample annotation table
# ============================================================

# Define output directory

output_directory <- file.path(
  "results",
  "tables"
)

# Create directory if it does not exist

if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

# Define output file

metadata_output <- file.path(
  output_directory,
  "integrated_sample_metadata.csv"
)

# Save integrated metadata

write.csv(
  integrated_metadata,
  file = metadata_output,
  row.names = FALSE
)

# Verify saved file

cat("\nSaved metadata file:\n")
print(metadata_output)

cat("\nFile size in bytes:\n")
print(file.info(metadata_output)$size)

# Verify number of saved records

saved_metadata <- read.csv(
  metadata_output,
  stringsAsFactors = FALSE
)

cat("\nSaved metadata dimensions:\n")
print(dim(saved_metadata))

# ============================================================
# STEP 19: Investigate unusually low-count sample
# ============================================================

# Define the sample requiring investigation

target_sample <- "IR9AD2L1"


# Check whether the sample exists

stopifnot(
  target_sample %in% colnames(counts)
)


# Calculate library sizes for all samples

library_sizes <- colSums(
  counts
)


# Calculate detected expression records

detected_records <- colSums(
  counts > 0
)


# Create sample-level QC table

sample_qc <- data.frame(
  
  library_name = names(library_sizes),
  
  library_size = as.numeric(library_sizes),
  
  detected_records = as.numeric(detected_records),
  
  stringsAsFactors = FALSE
  
)


# Display QC information for the target sample

cat("\nTarget sample quality information:\n")

print(
  sample_qc[
    sample_qc$library_name == target_sample,
  ]
)


# Display overall sample QC statistics

cat("\nLibrary size summary:\n")

print(
  summary(sample_qc$library_size)
)


cat("\nDetected expression records summary:\n")

print(
  summary(sample_qc$detected_records)
)

# ============================================================
# STEP 19.2: Compare sample with experimental group
# ============================================================


# Identify target sample metadata

target_metadata <- integrated_metadata[
  integrated_metadata$library_name == target_sample,
]


cat("\nTarget sample metadata:\n")

print(target_metadata)


# Identify samples from the same experimental group

same_group <- integrated_metadata[
  
  integrated_metadata$treatment == target_metadata$treatment &
    
    integrated_metadata$day == target_metadata$day &
    
    integrated_metadata$tissue == target_metadata$tissue,
  
]


cat("\nSamples from the same experimental group:\n")

print(
  same_group[
    ,
    c(
      "geo_accession",
      "library_name",
      "treatment",
      "day",
      "tissue",
      "replicate"
    )
  ]
)


# Extract QC metrics for the same group

group_qc <- merge(
  
  same_group,
  
  sample_qc,
  
  by = "library_name"
  
)


cat("\nQuality metrics for the experimental group:\n")

print(
  
  group_qc[
    ,
    c(
      "library_name",
      "replicate",
      "library_size",
      "detected_records"
    )
  ]
  
)

# ============================================================
# STEP 19.3: Calculate relative library size
# ============================================================

# Calculate median library size

median_library_size <- median(
  sample_qc$library_size
)


# Calculate target sample relative to median

target_library_size <- sample_qc$library_size[
  sample_qc$library_name == target_sample
]


relative_library_size <- (
  target_library_size / median_library_size
) * 100


cat("\nTarget sample library size:", target_library_size)

cat("\nMedian library size:", median_library_size)

cat(
  "\nTarget sample relative to median (%):",
  round(relative_library_size, 4),
  "%\n"
)


# Calculate fold difference

fold_difference <- median_library_size / target_library_size


cat(
  "\nMedian library size is",
  round(fold_difference, 2),
  "times larger than target sample.\n"
)

# ============================================================
# STEP 19.4: Visualize biological replicate comparison
# ============================================================

# Prepare plot data

group_plot_data <- group_qc[
  ,
  c(
    "library_name",
    "replicate",
    "library_size"
  )
]


# Create bar plot

barplot(
  height = group_plot_data$library_size,
  
  names.arg = group_plot_data$library_name,
  
  las = 2,
  
  ylab = "Total raw counts",
  
  main = "Library Size Comparison: Severe Drought Day 9",
  
  cex.names = 0.8
)


# Save figure

png(
  filename = file.path(
    "results",
    "figures",
    "severe_drought_day9_replicates.png"
  ),
  
  width = 1200,
  
  height = 800,
  
  res = 120
)

barplot(
  height = group_plot_data$library_size,
  
  names.arg = group_plot_data$library_name,
  
  las = 2,
  
  ylab = "Total raw counts",
  
  main = "Library Size Comparison: Severe Drought Day 9",
  
  cex.names = 0.8
)

dev.off()

# ============================================================
# STEP 19.5: Investigate target sample count distribution
# ============================================================

# Extract target sample counts

target_counts <- counts[[target_sample]]


# Summary statistics

cat("\nTarget sample count distribution:\n")

print(summary(target_counts))


# Count zero-expression records

zero_counts <- sum(target_counts == 0)

cat("\nZero-count records:", zero_counts)


# Count non-zero records

nonzero_counts <- sum(target_counts > 0)

cat("\nNon-zero records:", nonzero_counts)


# Identify highest expression records

top_counts <- sort(
  target_counts,
  decreasing = TRUE
)


cat("\nTop 10 expression counts:\n")

print(head(top_counts, 10))

# ============================================================
# STEP 19.6: Compare count distributions across replicates
# ============================================================

# Extract the three biological replicates

replicate_samples <- group_qc$library_name


# Calculate count distribution statistics

replicate_distribution <- data.frame(
  
  library_name = replicate_samples,
  
  total_counts = NA,
  
  zero_records = NA,
  
  nonzero_records = NA,
  
  maximum_count = NA,
  
  mean_count = NA
  
)


for (i in seq_along(replicate_samples)) {
  
  sample_name <- replicate_samples[i]
  
  current_counts <- counts[[sample_name]]
  
  replicate_distribution$total_counts[i] <- sum(current_counts)
  
  replicate_distribution$zero_records[i] <- sum(
    current_counts == 0
  )
  
  replicate_distribution$nonzero_records[i] <- sum(
    current_counts > 0
  )
  
  replicate_distribution$maximum_count[i] <- max(
    current_counts
  )
  
  replicate_distribution$mean_count[i] <- mean(
    current_counts
  )
  
}


# Display comparison

cat("\nBiological replicate count distribution:\n")

print(replicate_distribution)

# ============================================================
# STEP 19.7: Save replicate count distribution
# ============================================================

# Define output file

output_file <- file.path(
  "results",
  "tables",
  "replicate_count_distribution.csv"
)


# Save results

write.csv(
  replicate_distribution,
  file = output_file,
  row.names = FALSE
)


# Verify file creation

cat("\nFile successfully saved:",
    file.exists(output_file))


cat("\nOutput location:", output_file)