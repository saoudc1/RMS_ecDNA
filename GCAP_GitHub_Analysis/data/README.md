# Input data

Place the input files below in this directory before running the analysis.

## `Samples_for_validation.csv`

Required fields:

- `Sample_ID`
- `Reference_class`
- `GCAP_class`

Optional fields:

- `Reference_method` (`AA/WGS` or `FISH`)
- `Purity` or `Path_purity`

Accepted class labels include `ecDNA`/`circular`, `chromosomal`/`noncircular`, and `no focal`/`nofocal`.

## `RMS_GCAP_sample_info.csv`

Required fields:

- `sample`
- `class`
- `purity`

## `RMS_survival_with_purity.csv`

Required fields:

- `Sample.ID`
- `TumorPurity`

## `RMS_survival.csv`

Required fields:

- `Sample.ID`
- `Patient.ID`
- `Tumor_type`
- `Time_of_FU`
- a survival-status field such as `Status..NED..AWD.DOD.`

Disease-specific death is coded as `DOD` in the analysis.

## `gene_alterations.txt`

Tab-delimited file containing:

- `studyID:sampleId`
- `TP53`

## `ERMS_clinical.csv`

Required fields before R name normalization:

- `Patient ID`
- `Local Stage/Clinical Group`
- `age/sex`

The supplied analysis expects these to become `Patient.ID`, `Local.Stage.Clinical.Group`, and `age.sex` after `make.names()`.

