# CIN Signature Analysis in Rhabdomyosarcoma

This repository contains the R code used to quantify and compare copy-number instability (CIN) signatures in rhabdomyosarcoma using ASCAT-derived copy-number segments.

## Repository structure

```text
CIN_Signature_GitHub/
├── README.md
├── R/
│   └── CIN_signature_analysis.R
└── data/
    ├── ASCAT_segments/
    ├── Joined_final.csv
    ├── WGS_RMS_noMRN.csv
    └── README.md
```

## Requirements

The analysis uses the following R packages:

- CINSignatureQuantification
- dplyr
- tidyr
- stringr
- purrr
- tibble
- readr
- ggplot2
- ggsci

## Input data

Place the required files in the `data/` directory:

- `ASCAT_segments/`: ASCAT segment files matching `*_tumor.segments_raw.txt`; subdirectories are allowed.
- `Joined_final.csv`: AmpliconClassifier-derived sample and amplicon annotations.
- `WGS_RMS_noMRN.csv`: clinical/sample mapping file.

The ASCAT segment files must contain the columns `chr`, `startpos`, `endpos`, `nAraw`, and `nBraw`.

## Running the analysis

Run from the repository root:

```r
source("R/CIN_signature_analysis.R")
```

The script does not write files or save figures. Statistical result tables, processed datasets, and ggplot objects remain available in the R session.

Key result objects include:

- `activities`
- `merged`
- `merged_long`
- `ecdna_vs_non`
- `amp_global`
- `amp_pairwise`
- `subtype_global`
- `subtype_pairwise`
- `within_subtype_global`
- `within_subtype_pairwise`

Key figure objects include:

- `p_stack`
- `p_ecdna_binary`
- `p_amp`
- `p_subtype`

## Analysis overview

CIN signature activities are quantified from ASCAT allele-specific copy-number segments using `CINSignatureQuantification`. Samples are assigned to ecDNA, chromosomal amplification, or no focal amplification using a sample-level hierarchy of ecDNA > chromosomal amplification > no focal amplification. CIN signature exposures are compared across amplification classes and RMS subtypes using Wilcoxon rank-sum or Kruskal-Wallis tests as appropriate, with Benjamini-Hochberg correction for multiple testing.
