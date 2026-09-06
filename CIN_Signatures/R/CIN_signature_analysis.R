# =============================================================================
# ASCAT-derived CIN signature analysis in RMS
# Saoud et al. CCR 2026
# =============================================================================
suppressPackageStartupMessages({
  library(CINSignatureQuantification)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggsci)
})

# Inputs
base_dir <- file.path("data", "ASCAT_segments")
amplicon_file <- file.path("data", "Joined_final.csv")
clinical_file <- file.path("data", "WGS_RMS_noMRN.csv")

# =============================================================================
# 1. INPUTS
# =============================================================================


amp_levels <- c("no focal amplification", "chromosomal", "ecDNA")
amp_colors <- c(
  "no focal amplification" = "#008B45FF",
  "chromosomal" = "#EE0000FF",
  "ecDNA" = "#3B4992FF"
)

# =============================================================================
# 2. HELPERS
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

first_nonmissing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

as_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "1", "yes")
}

normalize_cin_sample <- function(x) {
  stringr::str_remove(as.character(x), stringr::regex("_[Tt]umor$"))
}

normalize_activities <- function(x) {
  x <- as.matrix(x)
  
  row_cx <- sum(stringr::str_detect(rownames(x) %||% character(), "^CX[0-9]+$"))
  col_cx <- sum(stringr::str_detect(colnames(x) %||% character(), "^CX[0-9]+$"))
  
  if (col_cx > row_cx) {
    out <- as.data.frame(x)
  } else if (row_cx > col_cx) {
    out <- as.data.frame(t(x))
  } else {
    stop("Could not determine orientation of CIN signature activities.")
  }
  
  tibble::rownames_to_column(out, "cin_sample_name")
}

safe_global_test <- function(dat, group_col) {
  
  dat <- dat |>
    dplyr::filter(
      !is.na(exposure),
      !is.na(.data[[group_col]])
    )
  
  groups <- unique(as.character(dat[[group_col]]))
  groups <- groups[!is.na(groups) & groups != ""]
  
  if (length(groups) < 2 || length(unique(dat$exposure)) < 2) {
    return(
      tibble::tibble(
        n = nrow(dat),
        n_groups = length(groups),
        test = NA_character_,
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  form <- stats::reformulate(group_col, response = "exposure")
  
  if (length(groups) == 2) {
    fit <- stats::wilcox.test(form, data = dat, exact = FALSE)
    test_name <- "Wilcoxon"
  } else {
    fit <- stats::kruskal.test(form, data = dat)
    test_name <- "Kruskal-Wallis"
  }
  
  tibble::tibble(
    n = nrow(dat),
    n_groups = length(groups),
    test = test_name,
    statistic = unname(fit$statistic),
    p_value = fit$p.value
  )
}

safe_pairwise_wilcox <- function(dat, group_col) {
  
  dat <- dat |>
    dplyr::filter(
      !is.na(exposure),
      !is.na(.data[[group_col]])
    )
  
  groups <- sort(unique(as.character(dat[[group_col]])))
  groups <- groups[!is.na(groups) & groups != ""]
  
  if (length(groups) < 2) return(tibble::tibble())
  
  pairs <- combn(groups, 2, simplify = FALSE)
  
  purrr::map_dfr(pairs, function(g) {
    
    x <- dat$exposure[as.character(dat[[group_col]]) == g[1]]
    y <- dat$exposure[as.character(dat[[group_col]]) == g[2]]
    
    if (length(x) == 0 || length(y) == 0 ||
        length(unique(c(x, y))) < 2) {
      
      stat <- NA_real_
      p <- NA_real_
      
    } else {
      
      fit <- stats::wilcox.test(x, y, exact = FALSE)
      stat <- unname(fit$statistic)
      p <- fit$p.value
    }
    
    tibble::tibble(
      group1 = g[1],
      group2 = g[2],
      n1 = length(x),
      n2 = length(y),
      statistic = stat,
      p_value = p
    )
  })
}

add_fdr_label <- function(df) {
  df |>
    dplyr::mutate(
      label = dplyr::case_when(
        is.na(FDR) ~ "",
        FDR <= 0.001 ~ "***",
        FDR <= 0.01 ~ "**",
        FDR <= 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

make_y_positions <- function(df) {
  pad <- max(df$exposure, na.rm = TRUE) * 0.05
  if (!is.finite(pad) || pad <= 0) pad <- 0.02
  
  df |>
    dplyr::group_by(signature) |>
    dplyr::summarise(
      y = max(exposure, na.rm = TRUE) + pad,
      .groups = "drop"
    )
}


# =============================================================================
# 3. READ ASCAT SEGMENTS
# =============================================================================

segment_files <- list.files(
  path = base_dir,
  pattern = "_tumor\\.segments_raw\\.txt$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(segment_files) == 0) {
  stop("No ASCAT *_tumor.segments_raw.txt files found in ", base_dir)
}

read_ascat_segments <- function(f) {
  
  x <- utils::read.delim(
    f,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required <- c(
    "chr",
    "startpos",
    "endpos",
    "nAraw",
    "nBraw"
  )
  
  missing <- setdiff(required, colnames(x))
  
  if (length(missing) > 0) {
    stop(
      "Missing required columns in ",
      basename(f),
      ": ",
      paste(missing, collapse = ", ")
    )
  }
  
  if (!"sample" %in% colnames(x)) {
    x$sample <- stringr::str_remove(
      basename(f),
      stringr::regex(
        "_tumor\\.segments_raw\\.txt$",
        ignore_case = TRUE
      )
    )
  }
  
  x
}

all_segments_df <- purrr::map_dfr(
  segment_files,
  read_ascat_segments
)

cin_input <- all_segments_df |>
  dplyr::mutate(
    chr = stringr::str_remove(
      as.character(chr),
      "^chr"
    ),
    startpos = as.numeric(startpos),
    endpos = as.numeric(endpos),
    nAraw = as.numeric(nAraw),
    nBraw = as.numeric(nBraw)
  ) |>
  dplyr::filter(
    chr %in% as.character(1:22),
    is.finite(startpos),
    is.finite(endpos),
    is.finite(nAraw),
    is.finite(nBraw),
    endpos >= startpos
  ) |>
  dplyr::transmute(
    chromosome = chr,
    start = startpos,
    end = endpos,
    segVal = nAraw + nBraw,
    sample = as.character(sample)
  )


# =============================================================================
# 4. QUANTIFY CIN SIGNATURES
# =============================================================================

result <- CINSignatureQuantification::quantifyCNSignatures(
  cin_input
)


activities_raw <- CINSignatureQuantification::getActivities(
  result
)

activities <- normalize_activities(
  activities_raw
)

sig_cols <- colnames(activities)[
  stringr::str_detect(
    colnames(activities),
    "^CX[0-9]+$"
  )
]

if (length(sig_cols) == 0) {
  stop("No CX signature columns detected.")
}

sig_cols <- sig_cols[
  order(
    as.numeric(
      stringr::str_extract(
        sig_cols,
        "[0-9]+"
      )
    )
  )
]

activities <- activities |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(sig_cols),
      as.numeric
    )
  )


# =============================================================================
# 5. OPTIONAL OVERVIEW STACKED BARPLOT
# =============================================================================

activity_mat <- as.matrix(
  activities[, sig_cols, drop = FALSE]
)

sample_summary <- tibble::tibble(
  cin_sample_name = activities$cin_sample_name,
  dominant_signature = sig_cols[
    max.col(
      activity_mat,
      ties.method = "first"
    )
  ],
  dominant_exposure = apply(
    activity_mat,
    1,
    max,
    na.rm = TRUE
  )
) |>
  dplyr::arrange(
    dominant_signature,
    dplyr::desc(dominant_exposure)
  )

sample_order <- sample_summary$cin_sample_name

act_long <- activities |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(sig_cols),
    names_to = "signature",
    values_to = "exposure"
  ) |>
  dplyr::mutate(
    signature = factor(
      signature,
      levels = sig_cols
    )
  )

p_stack <- ggplot2::ggplot(
  act_long,
  ggplot2::aes(
    x = factor(
      cin_sample_name,
      levels = sample_order
    ),
    y = exposure,
    fill = signature
  )
) +
  ggsci::scale_fill_igv() +
  ggplot2::geom_col(width = 0.9) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    x = NULL,
    y = "CIN signature exposure",
    fill = "Signature"
  )


# =============================================================================
# 6. BUILD SAMPLE METADATA
# =============================================================================

amp <- readr::read_csv(
  amplicon_file,
  show_col_types = FALSE
)

clinical <- utils::read.csv(
  clinical_file,
  stringsAsFactors = FALSE,
  check.names = TRUE
)

# Remove accidental leading/trailing whitespace or UTF-8 BOM from headers.
colnames(clinical) <- trimws(colnames(clinical))
colnames(clinical)[1] <- sub("^\\ufeff", "", colnames(clinical)[1])

required_amp <- c(
  "Sample.name",
  "Classification",
  "cancer_subtype",
  "patient"
)

required_clinical <- c(
  "DMP",
  "Sample.System.ID"
)

if (length(setdiff(required_amp, colnames(amp))) > 0) {
  stop(
    "Missing required amp columns: ",
    paste(
      setdiff(required_amp, colnames(amp)),
      collapse = ", "
    )
  )
}

if (length(setdiff(required_clinical, colnames(clinical))) > 0) {
  stop(
    "Missing required clinical columns: ",
    paste(
      setdiff(required_clinical, colnames(clinical)),
      collapse = ", "
    )
  )
}

# ---------------------------------------------------------------------------
# Sample-level amplification status directly from AmpliconClassifier output.
# A sample may have multiple feature rows, so collapse to one row per sample.
# ---------------------------------------------------------------------------

amp_summary <- amp |>
  dplyr::mutate(
    Classification = as.character(Classification)
  ) |>
  dplyr::group_by(
    Sample.name
  ) |>
  dplyr::summarise(
    patient_id = first_nonmissing(patient),
    cancer_subtype = first_nonmissing(cancer_subtype),
    
    has_ecDNA = any(
      Classification == "ecDNA",
      na.rm = TRUE
    ),
    
    has_chromosomal = any(
      Classification %in% c(
        "Linear",
        "Complex-non-cyclic",
        "BFB",
        "FAN"
      ),
      na.rm = TRUE
    ),
    
    has_no_focal = any(
      Classification == "No focal amplification",
      na.rm = TRUE
    ),
    
    amp_status = dplyr::case_when(
      has_ecDNA ~ "ecDNA",
      has_chromosomal ~ "chromosomal",
      has_no_focal ~ "no focal amplification",
      TRUE ~ NA_character_
    ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    amp_status = factor(
      amp_status,
      levels = amp_levels
    )
  )

unclassified_amp_samples <- amp_summary |>
  dplyr::filter(
    is.na(amp_status)
  )

if (nrow(unclassified_amp_samples) > 0) {
  
  warning(
    nrow(unclassified_amp_samples),
    " amp samples could not be assigned an amplification class."
  )
}

# ---------------------------------------------------------------------------
# clinical is used ONLY for Sample.name -> Sample.System.ID mapping.
# ---------------------------------------------------------------------------

clinical_map_check <- clinical |>
  dplyr::filter(
    !is.na(DMP),
    DMP != "",
    !is.na(Sample.System.ID),
    Sample.System.ID != ""
  ) |>
  dplyr::distinct(
    DMP,
    Sample.System.ID
  ) |>
  dplyr::count(
    DMP,
    name = "n_system_ids"
  )

bad_dmp_map <- clinical_map_check |>
  dplyr::filter(
    n_system_ids > 1
  )

if (nrow(bad_dmp_map) > 0) {
  
  stop(
    "Some DMP identifiers map to more than one Sample.System.ID."
  )
}

clinical_map <- clinical |>
  dplyr::filter(
    !is.na(DMP),
    DMP != "",
    !is.na(Sample.System.ID),
    Sample.System.ID != ""
  ) |>
  dplyr::group_by(
    DMP
  ) |>
  dplyr::summarise(
    Sample.System.ID = first_nonmissing(
      Sample.System.ID
    ),
    .groups = "drop"
  )

sample_metadata <- amp_summary |>
  dplyr::left_join(
    clinical_map,
    by = c(
      "Sample.name" = "DMP"
    )
  ) |>
  dplyr::select(
    Sample.name,
    Sample.System.ID,
    patient_id,
    cancer_subtype,
    amp_status,
    has_ecDNA,
    has_chromosomal,
    has_no_focal
  )

# ---------------------------------------------------------------------------
# Match ASCAT/CIN sample names to clinical Sample.System.ID.
# ---------------------------------------------------------------------------

cin_metadata <- activities |>
  dplyr::transmute(
    cin_sample_name,
    Sample.System.ID = normalize_cin_sample(
      cin_sample_name
    )
  ) |>
  dplyr::left_join(
    sample_metadata,
    by = "Sample.System.ID"
  )

unmatched <- cin_metadata |>
  dplyr::filter(
    is.na(amp_status)
  )


if (nrow(unmatched) > 0) {
  warning(
    nrow(unmatched),
    " CIN samples did not match amp through the clinical DMP mapping."
  )
}

merged <- cin_metadata |>
  dplyr::filter(
    !is.na(amp_status)
  ) |>
  dplyr::left_join(
    activities,
    by = "cin_sample_name"
  )

if (nrow(merged) == 0) {
  stop("No CIN samples matched amp/clinical metadata.")
}


# Sample counts used in the analyses.
sample_counts_amp <- merged |>
  dplyr::count(
    amp_status,
    name = "n_samples"
  )

sample_counts_subtype <- merged |>
  dplyr::count(
    cancer_subtype,
    name = "n_samples"
  )

sample_counts_cross <- merged |>
  dplyr::count(
    cancer_subtype,
    amp_status,
    name = "n_samples"
  )




# Repeated-patient sample counts
patient_counts <- merged |>
  dplyr::filter(
    !is.na(patient_id),
    patient_id != ""
  ) |>
  dplyr::distinct(
    patient_id,
    Sample.System.ID
  ) |>
  dplyr::count(
    patient_id,
    name = "n_samples"
  ) |>
  dplyr::arrange(
    dplyr::desc(n_samples)
  )


# Long-format dataset used for all signature comparisons.
merged_long <- merged |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(sig_cols),
    names_to = "signature",
    values_to = "exposure"
  ) |>
  dplyr::mutate(
    signature = factor(
      signature,
      levels = sig_cols
    ),
    exposure = as.numeric(
      exposure
    ),
    amp_status = factor(
      amp_status,
      levels = amp_levels
    )
  )

# =============================================================================
# 7. ecDNA VS NON-ecDNA
# =============================================================================

ecdna_binary_long <- merged_long |>
  dplyr::mutate(
    ecDNA_status = dplyr::if_else(
      amp_status == "ecDNA",
      "ecDNA",
      "non-ecDNA"
    ),
    ecDNA_status = factor(
      ecDNA_status,
      levels = c(
        "non-ecDNA",
        "ecDNA"
      )
    )
  )

ecdna_vs_non <- ecdna_binary_long |>
  dplyr::group_by(
    signature
  ) |>
  dplyr::group_modify(
    ~ safe_global_test(
      .x,
      "ecDNA_status"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::arrange(
    FDR,
    p_value
  )



ecdna_vs_non_anno <- ecdna_vs_non |>
  add_fdr_label() |>
  dplyr::left_join(
    make_y_positions(ecdna_binary_long),
    by = "signature"
  ) |>
  dplyr::mutate(
    signature = factor(signature, levels = sig_cols)
  )

p_ecdna_binary <- ggplot2::ggplot(
  ecdna_binary_long,
  ggplot2::aes(
    x = signature,
    y = exposure,
    fill = ecDNA_status
  )
) +
  ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = ecDNA_status),
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.2,
    alpha = 0.7
  ) +
  ggplot2::geom_text(
    data = ecdna_vs_non_anno,
    ggplot2::aes(
      x = signature,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "non-ecDNA" = "#BDBDBD",
      "ecDNA" = "#3B4992FF"
    )
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "non-ecDNA" = "#BDBDBD",
      "ecDNA" = "#3B4992FF"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.02, 0.12))
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "right"
  ) +
  ggplot2::labs(
    x = "CIN Signature",
    y = "Exposure",
    fill = "Amplicon status",
    color = "Amplicon status"
  )


# =============================================================================
# 8. ecDNA VS CHROMOSOMAL VS NO FOCAL AMPLIFICATION
# =============================================================================

amp_global <- merged_long |>
  dplyr::group_by(
    signature
  ) |>
  dplyr::group_modify(
    ~ safe_global_test(
      .x,
      "amp_status"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::arrange(
    FDR,
    p_value
  )

amp_pairwise <- merged_long |>
  dplyr::group_by(
    signature
  ) |>
  dplyr::group_modify(
    ~ safe_pairwise_wilcox(
      .x,
      "amp_status"
    )
  ) |>
  dplyr::group_by(
    signature
  ) |>
  dplyr::mutate(
    FDR_within_signature = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR_all_pairwise = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::left_join(
    amp_global |>
      dplyr::select(
        signature,
        global_FDR = FDR
      ),
    by = "signature"
  ) |>
  dplyr::arrange(
    global_FDR,
    FDR_within_signature
  )



amp_anno <- amp_global |>
  add_fdr_label() |>
  dplyr::left_join(
    make_y_positions(merged_long),
    by = "signature"
  ) |>
  dplyr::mutate(
    signature = factor(signature, levels = sig_cols)
  )

p_amp <- ggplot2::ggplot(
  merged_long,
  ggplot2::aes(
    x = signature,
    y = exposure,
    fill = amp_status
  )
) +
  ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = amp_status),
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.2,
    alpha = 0.7
  ) +
  ggplot2::geom_text(
    data = amp_anno,
    ggplot2::aes(
      x = signature,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(
    values = amp_colors,
    drop = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = amp_colors,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.02, 0.12))
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "right"
  ) +
  ggplot2::labs(
    x = "CIN Signature",
    y = "Exposure",
    fill = "Amplicon class",
    color = "Amplicon class"
  )


# =============================================================================
# 9. CIN SIGNATURES BY RMS SUBTYPE
# =============================================================================

subtype_long <- merged_long |>
  dplyr::filter(
    !is.na(cancer_subtype),
    cancer_subtype != ""
  )

subtype_global <- subtype_long |>
  dplyr::group_by(
    signature
  ) |>
  dplyr::group_modify(
    ~ safe_global_test(
      .x,
      "cancer_subtype"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::arrange(
    FDR,
    p_value
  )

subtype_pairwise <- subtype_long |>
  dplyr::group_by(
    signature
  ) |>
  dplyr::group_modify(
    ~ safe_pairwise_wilcox(
      .x,
      "cancer_subtype"
    )
  ) |>
  dplyr::group_by(
    signature
  ) |>
  dplyr::mutate(
    FDR_within_signature = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR_all_pairwise = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::left_join(
    subtype_global |>
      dplyr::select(
        signature,
        global_FDR = FDR
      ),
    by = "signature"
  )



subtype_anno <- subtype_global |>
  add_fdr_label() |>
  dplyr::left_join(
    make_y_positions(subtype_long),
    by = "signature"
  ) |>
  dplyr::mutate(
    signature = factor(signature, levels = sig_cols)
  )

p_subtype <- ggplot2::ggplot(
  subtype_long,
  ggplot2::aes(
    x = signature,
    y = exposure,
    fill = cancer_subtype
  )
) +
  ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = NA
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = cancer_subtype),
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.2,
    alpha = 0.7
  ) +
  ggplot2::geom_text(
    data = subtype_anno,
    ggplot2::aes(
      x = signature,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  ggsci::scale_fill_igv() +
  ggsci::scale_color_igv() +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.02, 0.12))
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "right"
  ) +
  ggplot2::labs(
    x = "CIN Signature",
    y = "Exposure",
    fill = "RMS subtype",
    color = "RMS subtype"
  )


# =============================================================================
# 10. AMPLICON-CLASS DIFFERENCES WITHIN EACH RMS SUBTYPE
# =============================================================================

within_subtype_global <- subtype_long |>
  dplyr::group_by(
    cancer_subtype,
    signature
  ) |>
  dplyr::group_modify(
    ~ safe_global_test(
      .x,
      "amp_status"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(
    cancer_subtype
  ) |>
  dplyr::mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(
    cancer_subtype,
    FDR,
    p_value
  )

within_subtype_pairwise <- subtype_long |>
  dplyr::group_by(
    cancer_subtype,
    signature
  ) |>
  dplyr::group_modify(
    ~ safe_pairwise_wilcox(
      .x,
      "amp_status"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(
    cancer_subtype,
    signature
  ) |>
  dplyr::mutate(
    FDR_within_signature = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(
    cancer_subtype
  ) |>
  dplyr::mutate(
    FDR_within_subtype = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    within_subtype_global |>
      dplyr::select(
        cancer_subtype,
        signature,
        global_FDR = FDR
      ),
    by = c(
      "cancer_subtype",
      "signature"
    )
  )



# One grouped boxplot per subtype, with all CIN signatures on the x-axis.
subtypes <- sort(
  unique(
    as.character(
      subtype_long$cancer_subtype
    )
  )
)

for (st in subtypes) {
  
  plot_data <- subtype_long |>
    dplyr::filter(
      cancer_subtype == st
    )
  
  if (dplyr::n_distinct(plot_data$amp_status) < 2) {
    next
  }
  
  plot_anno <- within_subtype_global |>
    dplyr::filter(
      cancer_subtype == st
    ) |>
    add_fdr_label() |>
    dplyr::left_join(
      make_y_positions(plot_data),
      by = "signature"
    ) |>
    dplyr::mutate(
      signature = factor(
        signature,
        levels = sig_cols
      )
    )
  
  p_within <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = signature,
      y = exposure,
      fill = amp_status
    )
  ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7,
      outlier.shape = NA
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = amp_status),
      position = ggplot2::position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.8
      ),
      size = 1.2,
      alpha = 0.7
    ) +
    ggplot2::geom_text(
      data = plot_anno,
      ggplot2::aes(
        x = signature,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.2
    ) +
    ggplot2::scale_fill_manual(
      values = amp_colors,
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = amp_colors,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.12))
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = st,
      x = "CIN Signature",
      y = "Exposure",
      fill = "Amplicon class",
      color = "Amplicon class"
    )
  
}

