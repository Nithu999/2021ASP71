# Stage 9: Inspect experimental groups

metadata <- read.csv(
  "results/tables/integrated_sample_metadata.csv"
)

cat("\nMetadata dimensions:\n")
print(dim(metadata))

cat("\nMetadata column names:\n")
print(names(metadata))

cat("\nMetadata structure:\n")
str(metadata)

cat("\nFirst six rows:\n")
print(head(metadata))

cat("\nExperimental group counts:\n")
print(table(metadata$experimental_group, useNA = "ifany"))


# STEP 9.1: Inspect experimental conditions

cat("\nTreatment groups:\n")
print(table(metadata$treatment))

cat("\nTime points:\n")
print(table(metadata$day))

cat("\nTissue types:\n")
print(table(metadata$tissue))

cat("\nTreatment and time-point combinations:\n")
print(with(metadata, table(treatment, day)))

cat("\nTreatment, time-point and tissue combinations:\n")
print(with(metadata, table(treatment, day, tissue)))

cat("\nReplicate distribution:\n")
print(with(metadata, table(treatment, day, replicate)))