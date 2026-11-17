# CropStress Intelligence — organized R workflow

## Folder layout

```text
CropStress_Intelligence/
├── code/
│   ├── 00_run_all.R
│   ├── 01_extract_metadata.R
│   ├── 02_import_counts.R
│   ├── 03_integrate_metadata.R
│   ├── 04_sample_quality_control.R
│   └── 05_gene_identifier_inspection.R
├── data/
│   └── raw/
├── references/
└── results/
    ├── figures/
    └── tables/
```

## How to use

1. Copy the files from this `code/` folder into your project's existing `code/` folder.
2. Keep the raw count file in `data/raw/`.
3. Open the RStudio project so the project root is the working directory.
4. Open and run `code/00_run_all.R`.

The scripts run in sequence and share objects in the same R environment. Run the main runner, not the individual stages out of order.

## Important notes

- This is a structural split of the existing `01_dataset_metadata.R` workflow. The original analysis blocks are preserved rather than redesigned.
- The sample flagged as `IR9AD2L1` is investigated but is not automatically removed.
- The current gene-identifier inspection checks R row names. If the count file has no gene-ID column, those row names may simply be row numbers. Do not treat them as gene identifiers until the original file format or a matching annotation file confirms this.
- The workflow still represents dataset preparation and initial quality control. The ICA3 also requires exploratory analysis, at least one suitable main analytical method, interpretation, and a full report. These stages must still be added.
- Before submission, check that all scripts run from start to finish on your machine and that the results are reproducible.
