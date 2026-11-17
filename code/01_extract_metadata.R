# CropStress Intelligence
# Stage 1: Download and extract GEO sample metadata
# Input: GSE288615 GEO SOFT family file
# Output: sample_metadata object in the R environment

# STEP 1: Define official GEO download URL

geo_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE288nnn/GSE288615/soft/",
  "GSE288615_family.soft.gz"
)


# STEP 2: Define metadata output location
metadata_file <- file.path(
  "references",
  "GSE288615_family.soft.gz"
)


# STEP 3: Check references directory

if (!dir.exists("references")) {
  dir.create("references")
}


# STEP 4: Download official GEO metadata

download.file(
  url = geo_url,
  destfile = metadata_file,
  mode = "wb"
)


# STEP 5: Verify downloaded file

cat("\nDownloaded metadata file size:\n")

print(file.info(metadata_file)$size)


# STEP 6: Read compressed GEO metadata

geo_connection <- gzfile(
  metadata_file,
  open = "rt"
)

soft_lines <- readLines(
  geo_connection,
  warn = FALSE
)

close(geo_connection)


# STEP 7: Count official GEO sample records

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


# STEP 8: Inspect sample annotation fields

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


# STEP 9: Extract individual GEO sample records


sample_starts <- grep(
  "^\\^SAMPLE =",
  soft_lines
)


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

      trimws(
        sub(
          paste0("^.*", key, ":\\s*"),
          "",
          result[1]
        )
      )
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


sample_metadata <- do.call(
  rbind,
  sample_metadata_list
)


cat("\nFirst six sample records:\n")

print(head(sample_metadata))


cat("\nMetadata dimensions:\n")

print(dim(sample_metadata))


# STEP 10: Validate extracted sample metadata

cat("\nMetadata validation:\n")


cat(
  "Total samples:",
  nrow(sample_metadata),
  "\n"
)


cat(
  "Missing values:",
  sum(is.na(sample_metadata)),
  "\n"
)


cat(
  "Duplicate GEO accessions:",
  sum(duplicated(sample_metadata$geo_accession)),
  "\n"
)


cat(
  "Unique GEO accessions:",
  length(unique(sample_metadata$geo_accession)),
  "\n"
)


# STEP 11: Summarise experimental design


cat("\nTreatment groups:\n")

print(unique(sample_metadata$treatment))


# Preserve original labels such as "0 (2h)"

cat("\nCollection time points:\n")

print(unique(sample_metadata$day))


cat("\nTissue types:\n")

print(unique(sample_metadata$tissue))


cat("\nSamples per treatment:\n")

print(table(sample_metadata$treatment))


cat("\nExperimental group distribution:\n")

group_distribution <- xtabs(
  ~ treatment + day + tissue,
  data = sample_metadata
)

print(group_distribution)


