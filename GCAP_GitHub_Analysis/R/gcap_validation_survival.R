# =============================================================================
# GCAP validation, purity sensitivity, and survival analyses
# =============================================================================
# Run this script from the repository root.
# Input files are expected in ./data/ (see README.md).

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(survival)
  library(coxphf)
})

# =============================================================================
# 1. INPUT FILES AND SETTINGS
# =============================================================================

data_dir <- "data"

validation_file    <- file.path(data_dir, "Samples_for_validation.csv")
gcap_file          <- file.path(data_dir, "RMS_GCAP_sample_info.csv")
purity_file        <- file.path(data_dir, "RMS_survival_with_purity.csv")
survival_file      <- file.path(data_dir, "RMS_survival.csv")
gene_file          <- file.path(data_dir, "gene_alterations.txt")
erms_clinical_file <- file.path(data_dir, "ERMS_clinical.csv")

cutoffs <- c(NA_real_, 20, 30, 40, 50)
cutoff_levels <- c("No cutoff", "20%", "30%", "40%", "50%")

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

cutoff_label <- function(x) {
  ifelse(is.na(x), "No cutoff", paste0(x, "%"))
}

first_nonmissing_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else x[1]
}

first_nonmissing_character <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[1]
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

safe_divide <- function(x, y) {
  if (is.na(y) || y == 0) NA_real_ else x / y
}

exact_ci <- function(success, total) {
  if (is.na(total) || total == 0) return(c(NA_real_, NA_real_))
  unname(stats::binom.test(success, total)$conf.int)
}

clean_class <- function(x) {
  x <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  x <- stringr::str_replace_all(x, "_", " ")
  x <- stringr::str_replace_all(x, "-", " ")
  
  dplyr::case_when(
    x %in% c("ecdna", "ec dna", "circular") ~ "ecDNA",
    x %in% c("chromosomal", "noncircular", "non circular") ~ "chromosomal",
    x %in% c("nofocal", "no focal", "no amplification") ~ "no focal",
    is.na(x) | x == "" ~ NA_character_,
    TRUE ~ x
  )
}

as_binary_altered <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))
  x_num <- suppressWarnings(as.numeric(x_chr))
  out <- rep(NA_integer_, length(x_chr))
  
  out[!is.na(x_num)] <- as.integer(x_num[!is.na(x_num)] != 0)
  
  out[is.na(out) & x_chr %in% c(
    "true", "t", "yes", "y", "altered", "mutated", "mutation",
    "positive", "pos", "amp", "amplified"
  )] <- 1L
  
  out[is.na(out) & x_chr %in% c(
    "false", "f", "no", "n", "wt", "wild type", "wildtype",
    "negative", "neg", ""
  )] <- 0L
  
  out
}

# =============================================================================
# 3. LOAD DATA
# =============================================================================

gcap_results <- read.csv(
  gcap_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

path_purity <- read.csv(
  purity_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# =============================================================================
# 4. GCAP VALIDATION COHORT
# =============================================================================

validation_raw <- read.csv(
  validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# The validation file is expected to contain Sample_ID, Reference_class,
# and GCAP_class. 
if (!all(c("Sample_ID", "Reference_class", "GCAP_class") %in% names(validation_raw))) {
  if (ncol(validation_raw) < 4) {
    stop("Validation file must contain Sample_ID, Reference_class, and GCAP_class.")
  }
  colnames(validation_raw)[2:4] <- c("Sample_ID", "Reference_class", "GCAP_class")
}

purity_cols <- which(
  grepl(
    "^(Purity|Path_purity|Path\\.purity|Path purity)$",
    trimws(names(validation_raw)),
    ignore.case = TRUE
  )
)

if (length(purity_cols) > 0) {
  sheet_purity_matrix <- validation_raw[, purity_cols, drop = FALSE]
  Sheet_purity_raw <- apply(
    sheet_purity_matrix,
    1,
    function(x) {
      z <- suppressWarnings(as.numeric(trimws(as.character(x))))
      z <- z[is.finite(z)]
      if (length(z) == 0) NA_real_ else z[1]
    }
  )
} else {
  Sheet_purity_raw <- rep(NA_real_, nrow(validation_raw))
}

Sheet_purity_percent <- dplyr::case_when(
  is.na(Sheet_purity_raw) ~ NA_real_,
  Sheet_purity_raw <= 1 ~ Sheet_purity_raw * 100,
  TRUE ~ Sheet_purity_raw
)

# Reference method can be supplied directly. Otherwise FISH cases are detected
# from the row text, matching the original analysis logic.
if ("Reference_method" %in% names(validation_raw)) {
  Reference_method <- as.character(validation_raw$Reference_method)
} else {
  FISH_from_sheet <- apply(
    validation_raw,
    1,
    function(x) any(grepl("FISH", as.character(x), ignore.case = TRUE), na.rm = TRUE)
  )
  Reference_method <- ifelse(FISH_from_sheet, "FISH", "AA/WGS")
}

validation <- tibble::tibble(
  Sample_ID = trimws(as.character(validation_raw$Sample_ID)),
  Reference_class = clean_class(validation_raw$Reference_class),
  GCAP_class = clean_class(validation_raw$GCAP_class),
  Reference_method = Reference_method,
  Sheet_purity_percent = Sheet_purity_percent
) |>
  dplyr::mutate(
    Reference_ecDNA = ifelse(Reference_class == "ecDNA", "Positive", "Negative"),
    GCAP_ecDNA = ifelse(GCAP_class == "ecDNA", "Positive", "Negative")
  )

path_purity_map <- path_purity |>
  dplyr::transmute(
    Sample_ID = trimws(as.character(Sample.ID)),
    File_purity_percent = suppressWarnings(as.numeric(TumorPurity))
  ) |>
  dplyr::group_by(Sample_ID) |>
  dplyr::summarise(
    File_purity_percent = first_nonmissing_numeric(File_purity_percent),
    .groups = "drop"
  )

gcap_purity_map <- gcap_results |>
  dplyr::transmute(
    Sample_ID = trimws(as.character(sample)),
    GCAP_purity_raw = suppressWarnings(as.numeric(purity))
  ) |>
  dplyr::group_by(Sample_ID) |>
  dplyr::summarise(
    GCAP_purity_raw = first_nonmissing_numeric(GCAP_purity_raw),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    GCAP_purity_percent = dplyr::case_when(
      is.na(GCAP_purity_raw) ~ NA_real_,
      GCAP_purity_raw <= 1 ~ GCAP_purity_raw * 100,
      TRUE ~ GCAP_purity_raw
    )
  )

validation_purity <- validation |>
  dplyr::left_join(path_purity_map, by = "Sample_ID") |>
  dplyr::left_join(gcap_purity_map, by = "Sample_ID") |>
  dplyr::mutate(
    Path_purity_percent = dplyr::coalesce(
      File_purity_percent,
      Sheet_purity_percent
    ),
    # For FISH-reference cases, use GCAP-estimated purity when available.
    Path_purity_percent = dplyr::case_when(
      Reference_method == "FISH" & !is.na(GCAP_purity_percent) ~ GCAP_purity_percent,
      TRUE ~ Path_purity_percent
    )
  )

calc_performance <- function(df, cutoff) {
  TP <- sum(df$Reference_ecDNA == "Positive" & df$GCAP_ecDNA == "Positive", na.rm = TRUE)
  FN <- sum(df$Reference_ecDNA == "Positive" & df$GCAP_ecDNA == "Negative", na.rm = TRUE)
  FP <- sum(df$Reference_ecDNA == "Negative" & df$GCAP_ecDNA == "Positive", na.rm = TRUE)
  TN <- sum(df$Reference_ecDNA == "Negative" & df$GCAP_ecDNA == "Negative", na.rm = TRUE)
  
  sens_ci <- exact_ci(TP, TP + FN)
  spec_ci <- exact_ci(TN, TN + FP)
  acc_ci  <- exact_ci(TP + TN, TP + TN + FP + FN)
  ppv_ci  <- exact_ci(TP, TP + FP)
  npv_ci  <- exact_ci(TN, TN + FN)
  
  tibble::tibble(
    cutoff = cutoff,
    cutoff_label = cutoff_label(cutoff),
    N = nrow(df),
    excluded = nrow(validation_purity) - nrow(df),
    TP = TP,
    FN = FN,
    FP = FP,
    TN = TN,
    Sensitivity = safe_divide(TP, TP + FN),
    Sensitivity_lower = sens_ci[1],
    Sensitivity_upper = sens_ci[2],
    Specificity = safe_divide(TN, TN + FP),
    Specificity_lower = spec_ci[1],
    Specificity_upper = spec_ci[2],
    Accuracy = safe_divide(TP + TN, TP + TN + FP + FN),
    Accuracy_lower = acc_ci[1],
    Accuracy_upper = acc_ci[2],
    PPV = safe_divide(TP, TP + FP),
    PPV_lower = ppv_ci[1],
    PPV_upper = ppv_ci[2],
    NPV = safe_divide(TN, TN + FN),
    NPV_lower = npv_ci[1],
    NPV_upper = npv_ci[2]
  )
}

validation_subset <- function(cutoff) {
  # Validation sensitivity analysis: the threshold is applied to all samples.
  if (is.na(cutoff)) return(validation_purity)
  
  validation_purity |>
    dplyr::filter(
      !is.na(Path_purity_percent),
      Path_purity_percent >= cutoff
    )
}

performance_results <- purrr::map_dfr(
  cutoffs,
  function(cutoff) calc_performance(validation_subset(cutoff), cutoff)
) |>
  dplyr::mutate(
    cutoff_label = factor(cutoff_label, levels = cutoff_levels)
  ) |>
  dplyr::arrange(cutoff_label)

validation_30_vs_50 <- performance_results |>
  dplyr::filter(!is.na(cutoff), cutoff %in% c(30, 50))

get_discordant_cases <- function(cutoff) {
  validation_subset(cutoff) |>
    dplyr::filter(
      !is.na(Reference_ecDNA),
      !is.na(GCAP_ecDNA),
      Reference_ecDNA != GCAP_ecDNA
    ) |>
    dplyr::mutate(
      cutoff = cutoff,
      cutoff_label = cutoff_label(cutoff),
      Discordance = dplyr::case_when(
        Reference_ecDNA == "Positive" & GCAP_ecDNA == "Negative" ~ "False negative",
        Reference_ecDNA == "Negative" & GCAP_ecDNA == "Positive" ~ "False positive"
      )
    ) |>
    dplyr::select(
      cutoff_label,
      Sample_ID,
      Discordance,
      Reference_method,
      Reference_class,
      GCAP_class,
      Path_purity_percent,
      GCAP_purity_percent
    )
}

discordant_no_cutoff <- get_discordant_cases(NA_real_)
discordant_30 <- get_discordant_cases(30)
discordant_50 <- get_discordant_cases(50)

# =============================================================================
# 5. SURVIVAL DATA
# =============================================================================

rms_survival_all <- read.csv(
  survival_file,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

status_candidates <- c(
  "Status..NED..AWD.DOD.",
  "Status",
  "status",
  "Survival.Status",
  "Vital.Status"
)

status_col <- status_candidates[status_candidates %in% names(rms_survival_all)][1]
if (is.na(status_col)) {
  stop("Could not identify the survival-status column.")
}

clinical_gcap <- rms_survival_all |>
  dplyr::transmute(
    Sample.ID = trimws(as.character(Sample.ID)),
    Patient.ID = trimws(as.character(Patient.ID)),
    Tumor_type = trimws(as.character(Tumor_type)),
    Status = trimws(as.character(.data[[status_col]])),
    Time_of_FU = suppressWarnings(as.numeric(Time_of_FU))
  ) |>
  dplyr::distinct()

path_purity_survival <- path_purity |>
  dplyr::transmute(
    sample = trimws(as.character(Sample.ID)),
    path_purity = suppressWarnings(as.numeric(TumorPurity))
  ) |>
  dplyr::group_by(sample) |>
  dplyr::summarise(
    path_purity = first_nonmissing_numeric(path_purity),
    .groups = "drop"
  )

gcap_survival_base <- gcap_results |>
  dplyr::transmute(
    sample = trimws(as.character(sample)),
    class = tolower(trimws(as.character(class))),
    computational_purity = suppressWarnings(as.numeric(purity))
  ) |>
  dplyr::left_join(path_purity_survival, by = "sample")

make_patient_survival <- function(cutoff) {
  if (is.na(cutoff)) {
    filtered <- gcap_survival_base
  } else {
    cutoff_fraction <- cutoff / 100
    
    filtered <- gcap_survival_base |>
      dplyr::mutate(
        passes_purity =
          !is.na(computational_purity) &
          computational_purity >= cutoff_fraction &
          !is.na(path_purity) &
          path_purity >= cutoff
      ) |>
      # Primary downstream rule:
      # retain circular/ecDNA calls regardless of purity; apply the threshold
      # to noncircular and no-focal calls.
      dplyr::filter(class == "circular" | passes_purity)
  }
  
  patient <- filtered |>
    dplyr::left_join(clinical_gcap, by = c("sample" = "Sample.ID")) |>
    dplyr::filter(
      !is.na(Patient.ID),
      class %in% c("circular", "noncircular", "nofocal")
    ) |>
    dplyr::group_by(Patient.ID) |>
    dplyr::summarise(
      class = dplyr::case_when(
        any(class == "circular", na.rm = TRUE) ~ "circular",
        any(class == "noncircular", na.rm = TRUE) ~ "noncircular",
        any(class == "nofocal", na.rm = TRUE) ~ "nofocal",
        TRUE ~ NA_character_
      ),
      Tumor_type = first_nonmissing_character(Tumor_type),
      Status = first_nonmissing_character(Status),
      Time_of_FU = safe_max(Time_of_FU),
      .groups = "drop"
    ) |>
    dplyr::filter(
      !is.na(class),
      !is.na(Status),
      !is.na(Time_of_FU),
      Time_of_FU >= 0
    ) |>
    dplyr::mutate(
      event = as.integer(Status == "DOD"),
      ecDNA_binary = factor(
        ifelse(class == "circular", "ecDNA", "non-ecDNA"),
        levels = c("non-ecDNA", "ecDNA")
      )
    )
  
  patient
}

run_overall_firth <- function(patient_data, cutoff) {
  df <- patient_data |>
    dplyr::filter(!is.na(ecDNA_binary)) |>
    droplevels()
  
  fit <- coxphf::coxphf(
    survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
    data = df
  )
  
  tibble::tibble(
    cutoff = cutoff,
    cutoff_label = cutoff_label(cutoff),
    n = nrow(df),
    deaths = sum(df$event, na.rm = TRUE),
    non_ecDNA_n = sum(df$ecDNA_binary == "non-ecDNA"),
    ecDNA_n = sum(df$ecDNA_binary == "ecDNA"),
    HR = exp(fit$coefficients[1]),
    CI_lower = fit$ci.lower[1],
    CI_upper = fit$ci.upper[1],
    P_value = fit$prob[1]
  )
}

# =============================================================================
# 6. ERMS MULTIVARIABLE ANALYSIS
# =============================================================================

genes_altered <- read.delim(
  gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_sample <- as.character(genes_altered[["studyID:sampleId"]])
gene_sample <- ifelse(grepl(":", gene_sample), sub("^.*:", "", gene_sample), gene_sample)

patient_tp53 <- tibble::tibble(
  Patient.ID = stringr::str_extract(gene_sample, "P-[0-9]+"),
  TP53 = as_binary_altered(genes_altered$TP53)
) |>
  dplyr::filter(!is.na(Patient.ID)) |>
  dplyr::group_by(Patient.ID) |>
  dplyr::summarise(
    TP53 = if (all(is.na(TP53))) NA_integer_ else as.integer(any(TP53 == 1, na.rm = TRUE)),
    .groups = "drop"
  )

erms_data <- read.csv(
  erms_clinical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(erms_data) <- make.names(names(erms_data), unique = TRUE)

stage_group_split <- stringr::str_split_fixed(
  trimws(as.character(erms_data$Local.Stage.Clinical.Group)),
  "/",
  2
)

age_sex_split <- stringr::str_split_fixed(
  trimws(as.character(erms_data$age.sex)),
  "/",
  2
)

erms_covariates <- tibble::tibble(
  Patient.ID = trimws(as.character(erms_data$Patient.ID)),
  Local_stage = suppressWarnings(as.numeric(trimws(stage_group_split[, 1]))),
  Clinical_group = trimws(stage_group_split[, 2]),
  Age = suppressWarnings(as.numeric(trimws(age_sex_split[, 1]))),
  Sex_code = toupper(trimws(age_sex_split[, 2]))
) |>
  dplyr::group_by(Patient.ID) |>
  dplyr::summarise(
    Local_stage = first_nonmissing_numeric(Local_stage),
    Clinical_group = first_nonmissing_character(Clinical_group),
    Age = first_nonmissing_numeric(Age),
    Sex_code = first_nonmissing_character(Sex_code),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Age_group = factor(
      dplyr::case_when(
        Age < 10 ~ "<10 years",
        Age >= 10 ~ ">=10 years",
        TRUE ~ NA_character_
      ),
      levels = c("<10 years", ">=10 years")
    ),
    Sex = factor(
      Sex_code,
      levels = c("F", "M"),
      labels = c("Female", "Male")
    ),
    Stage_group = factor(
      dplyr::case_when(
        Local_stage %in% c(1, 2) ~ "Stage 1-2",
        Local_stage %in% c(3, 4) ~ "Stage 3-4",
        TRUE ~ NA_character_
      ),
      levels = c("Stage 1-2", "Stage 3-4")
    ),
    Clinical_group_binary = factor(
      dplyr::case_when(
        Clinical_group %in% c("I", "II") ~ "Group I-II",
        Clinical_group %in% c("III", "IV") ~ "Group III-IV",
        TRUE ~ NA_character_
      ),
      levels = c("Group I-II", "Group III-IV")
    )
  )

run_erms_multivariable <- function(patient_data, cutoff) {
  df <- patient_data |>
    dplyr::filter(Tumor_type == "ERMS") |>
    dplyr::left_join(patient_tp53, by = "Patient.ID") |>
    dplyr::left_join(erms_covariates, by = "Patient.ID") |>
    dplyr::mutate(
      TP53_status = factor(
        TP53,
        levels = c(0, 1),
        labels = c("TP53 WT", "TP53 altered")
      )
    ) |>
    dplyr::filter(
      !is.na(ecDNA_binary),
      !is.na(TP53_status),
      !is.na(Age_group),
      !is.na(Sex),
      !is.na(Stage_group),
      !is.na(Clinical_group_binary),
      !is.na(Time_of_FU),
      !is.na(event)
    ) |>
    droplevels()
  
  fit <- coxphf::coxphf(
    survival::Surv(Time_of_FU, event) ~
      ecDNA_binary +
      TP53_status +
      Age_group +
      Sex +
      Stage_group +
      Clinical_group_binary,
    data = df,
    maxit = 1000,
    maxstep = 0.5
  )
  
  ec_idx <- grep("^ecDNA_binary", names(fit$coefficients))[1]
  
  tibble::tibble(
    cutoff = cutoff,
    cutoff_label = cutoff_label(cutoff),
    n = nrow(df),
    deaths = sum(df$event, na.rm = TRUE),
    HR = exp(fit$coefficients[ec_idx]),
    CI_lower = fit$ci.lower[ec_idx],
    CI_upper = fit$ci.upper[ec_idx],
    P_value = fit$prob[ec_idx]
  )
}

# =============================================================================
# 7. RUN PURITY SENSITIVITY ANALYSES
# =============================================================================

patient_survival_by_cutoff <- purrr::set_names(
  purrr::map(cutoffs, make_patient_survival),
  cutoff_levels
)

overall_results <- purrr::map2_dfr(
  patient_survival_by_cutoff,
  cutoffs,
  run_overall_firth
) |>
  dplyr::mutate(cutoff_label = factor(cutoff_label, levels = cutoff_levels)) |>
  dplyr::arrange(cutoff_label)

erms_results <- purrr::map2_dfr(
  patient_survival_by_cutoff,
  cutoffs,
  run_erms_multivariable
) |>
  dplyr::mutate(cutoff_label = factor(cutoff_label, levels = cutoff_levels)) |>
  dplyr::arrange(cutoff_label)

survival_retention <- purrr::map2_dfr(
  patient_survival_by_cutoff,
  cutoffs,
  function(df, cutoff) {
    tibble::tibble(
      cutoff = cutoff,
      cutoff_label = cutoff_label(cutoff),
      patients_retained = nrow(df),
      deaths = sum(df$event, na.rm = TRUE)
    )
  }
) |>
  dplyr::mutate(
    cutoff_label = factor(cutoff_label, levels = cutoff_levels),
    total_patients = patients_retained[is.na(cutoff)][1],
    patients_excluded = total_patients - patients_retained
  ) |>
  dplyr::arrange(cutoff_label)

survival_30_vs_50 <- survival_retention |>
  dplyr::filter(!is.na(cutoff), cutoff %in% c(30, 50))

