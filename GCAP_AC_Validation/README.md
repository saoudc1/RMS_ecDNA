# GCAP validation and survival analysis

This repository contains the R code used for GCAP validation, tumor-purity sensitivity analyses, and ecDNA survival analyses for Saoud et al RMS ecDNA manuscript.

## Repository structure

```text
GCAP_GitHub_Analysis/
├── README.md
├── .gitignore
├── R/
│   └── gcap_validation_survival.R
└── data/
    └── README.md
```

## Requirements

R packages:

```r
install.packages(c(
  "dplyr",
  "purrr",
  "stringr",
  "tibble",
  "survival",
  "coxphf"
))
```

## Running the analysis

Run R from the repository root and execute:

```r
source("R/gcap_validation_survival.R")
```

No output files are created. Results are retained in the R workspace.

## Main analysis objects

- `validation_purity`: validation cohort with tumor-purity information.
- `performance_results`: GCAP sensitivity, specificity, accuracy, PPV, and NPV across purity thresholds.
- `validation_30_vs_50`: direct comparison of the 30% and 50% validation thresholds.
- `discordant_no_cutoff`, `discordant_30`, `discordant_50`: discordant GCAP/reference cases.
- `patient_survival_by_cutoff`: patient-level survival datasets for each purity threshold.
- `survival_retention`: number of patients retained across purity thresholds.
- `overall_results`: overall RMS Firth Cox results for ecDNA versus non-ecDNA.
- `erms_results`: ERMS multivariable Firth Cox results across purity thresholds.
- `primary_validation_30`, `primary_overall_30`, `primary_erms_30`: results at the prespecified 30% threshold.

## Analysis definitions

For validation sensitivity analyses, the purity threshold is applied to all validation samples.

For downstream survival analyses, circular/ecDNA-positive calls are retained regardless of purity, whereas noncircular and no-focal calls are required to meet both the GCAP computational-purity and pathology-purity thresholds.

At the patient level, GCAP classification follows the hierarchy:

```text
ecDNA > chromosomal amplification > no focal amplification
```

The ERMS multivariable Firth Cox model adjusts for ecDNA status, TP53 alteration status, age group, sex, local stage, and clinical group, matching the supplied analysis script.

## Data availability

The repository does not include patient-level clinical or sequencing data by default. Place only data that are approved for public release in the `data/` directory. Do not upload protected health information or restricted-access genomic data to GitHub.

See `data/README.md` for the expected filenames and columns.
