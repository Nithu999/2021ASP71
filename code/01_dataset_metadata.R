
# ============================================================
# CropStress Intelligence
# Script: 01_dataset_metadata.R
# Purpose: Download official GEO metadata
# Dataset: GSE288615
# ============================================================

# Step 1: Define GEO download URL

geo_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE288nnn/GSE288615/soft/",
  "GSE288615_family.soft.gz"
)

# Step 2: Define output file path

metadata_file <- file.path(
  "references",
  "GSE288615_family.soft.gz"
)

# Step 3: Check references directory

if (!dir.exists("references")) {
  dir.create("references")
}

# Step 4: Download GEO metadata

download.file(
  url = geo_url,
  destfile = metadata_file,
  mode = "wb"
)

# Step 5: Verify downloaded file

file.info(metadata_file)$size


# Step 6: Read the compressed GEO metadata

geo_connection <- gzfile(metadata_file, open = "rt")

soft_lines <- readLines(
  geo_connection,
  warn = FALSE
)

close(geo_connection)

# Step 7: Count official GEO sample records

sample_records <- grep(
  "^\\^SAMPLE =",
  soft_lines,
  value = TRUE
)

cat(
  "Number of GEO sample records:",
  length(sample_records),
  "\n"
)

# Step 8: Inspect important sample annotation fields

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
# Step 9: Extract individual GEO sample records
# ============================================================

# Find the beginning of every sample record
sample_starts <- grep(
  "^\\^SAMPLE =",
  soft_lines
)

# Extract one sample's metadata
extract_field <- function(lines, pattern) {
  
  result <- grep(
    pattern,
    lines,
    value = TRUE
  )
  
  if (length(result) == 0) {
    return(NA_character_)
  }
  
  sub("^[^=]*=\\s*", "", result[1])
}

# Process all 72 sample records
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
    
    # Extract characteristics
    characteristics <- grep(
      "^!Sample_characteristics_ch1",
      sample_block,
      value = TRUE
    )
    
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
      
      trimws(sub(
        paste0("^.*", key, ":\\s*"),
        "",
        result[1]
      ))
    }
    
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

# Combine individual records into one table
sample_metadata <- do.call(
  rbind,
  sample_metadata_list
)

# Display the first records
head(sample_metadata)

# Check dimensions
dim(sample_metadata)


# ============================================================
# Step 10: Summarise the experimental design
# ============================================================

# 1. Identify all treatment groups

cat("\nTreatment groups:\n")
print(unique(sample_metadata$treatment))

# 2. Identify collection time points

cat("\nCollection days:\n")
print(sort(unique(as.numeric(sample_metadata$day))))

# 3. Identify tissue types

cat("\nTissue types:\n")
print(unique(sample_metadata$tissue))

# 4. Count samples by treatment

cat("\nSamples per treatment:\n")
print(table(sample_metadata$treatment))

# 5. Examine treatment, day and tissue combinations

cat("\nExperimental group distribution:\n")

group_distribution <- xtabs(
  ~ treatment + day + tissue,
  data = sample_metadata
)

print(group_distribution)
