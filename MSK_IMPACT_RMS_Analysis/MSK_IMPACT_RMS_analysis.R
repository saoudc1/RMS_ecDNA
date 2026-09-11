# =============================================================================
# MSK-IMPACT RMS Cohort
# Saoud et al.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

# 1. INPUT FILES
# Data are not included in this repository.
# Place the required input files in the R working directory or update the
# filenames below before running the analysis.

survival_file <- "RMS_ALL_survival_4_2026.csv"
gene_file <- "MDM2_CDK4_MYC_status.txt"
purity_file <- "RMS_ALL_survival_with_purity_4_2026.csv"
gcap_file <- "RMS_ALL_6_2026_sample_info.csv"
erms_clinical_file <- "ERMS_IMPACT_survival_data_LHW_no_MRN.csv"
recurrent_genes_file <- "RMS_6_2024_annotated.csv"
impact_gene_list_file <- "IMPACT_genes.csv"
fold_change_file <- "genes_fc.csv"

arms_oncoprint_file <- "ARMS_oncoprint.csv"
erms_oncoprint_file <- "ERMS_oncoprint_4_2026.csv"
ssrms_oncoprint_file <- "SSRMS_oncoprint_4_2026.csv"

gene_cols <- c("MDM2", "CDK4", "GLI1", "MYCN", "MYC", "TP53")
genes_12q <- c("MDM2", "CDK4", "GLI1")

tumor_order <- c("ARMS", "ERMS", "SSRMS")

tumor_colors <- c(
  "ARMS" = "#3F4E9A",
  "ERMS" = "#F00000",
  "SSRMS" = "#078C45"
)

# 2. HELPERS

run_fisher_one_vs_rest <- function(data, gene, tumor) {
  
  df <- data |>
    dplyr::mutate(
      gene_altered = .data[[gene]] == 1,
      in_tumor_type = Tumor_type == tumor
    )
  
  tab <- table(
    df$in_tumor_type,
    df$gene_altered
  )
  
  altered_in_tumor <- sum(
    df$gene_altered & df$in_tumor_type,
    na.rm = TRUE
  )
  
  total_in_tumor <- sum(
    df$in_tumor_type,
    na.rm = TRUE
  )
  
  altered_in_other <- sum(
    df$gene_altered & !df$in_tumor_type,
    na.rm = TRUE
  )
  
  total_in_other <- sum(
    !df$in_tumor_type,
    na.rm = TRUE
  )
  
  if (!all(dim(tab) == c(2, 2))) {
    return(
      tibble::tibble(
        Gene = gene,
        Tumor_type = tumor,
        altered_in_tumor = altered_in_tumor,
        total_in_tumor = total_in_tumor,
        altered_in_other = altered_in_other,
        total_in_other = total_in_other,
        percent_in_tumor = 100 * altered_in_tumor / total_in_tumor,
        percent_in_other = 100 * altered_in_other / total_in_other,
        odds_ratio = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- stats::fisher.test(tab)
  
  tibble::tibble(
    Gene = gene,
    Tumor_type = tumor,
    altered_in_tumor = altered_in_tumor,
    total_in_tumor = total_in_tumor,
    altered_in_other = altered_in_other,
    total_in_other = total_in_other,
    percent_in_tumor = 100 * altered_in_tumor / total_in_tumor,
    percent_in_other = 100 * altered_in_other / total_in_other,
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  )
}

run_fisher_12q <- function(data, tumor) {
  
  df <- data |>
    dplyr::mutate(
      in_tumor_type = Tumor_type == tumor,
      amp_12q = Amp_12q == 1
    )
  
  tab <- table(
    df$in_tumor_type,
    df$amp_12q
  )
  
  if (!all(dim(tab) == c(2, 2))) {
    return(
      tibble::tibble(
        Tumor_type = tumor,
        amplified_in_tumor = sum(df$amp_12q & df$in_tumor_type, na.rm = TRUE),
        total_in_tumor = sum(df$in_tumor_type, na.rm = TRUE),
        amplified_in_other = sum(df$amp_12q & !df$in_tumor_type, na.rm = TRUE),
        total_in_other = sum(!df$in_tumor_type, na.rm = TRUE),
        percent_in_tumor = NA_real_,
        percent_in_other = NA_real_,
        odds_ratio = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  test <- stats::fisher.test(tab)
  
  amplified_in_tumor <- sum(
    df$amp_12q & df$in_tumor_type,
    na.rm = TRUE
  )
  
  total_in_tumor <- sum(
    df$in_tumor_type,
    na.rm = TRUE
  )
  
  amplified_in_other <- sum(
    df$amp_12q & !df$in_tumor_type,
    na.rm = TRUE
  )
  
  total_in_other <- sum(
    !df$in_tumor_type,
    na.rm = TRUE
  )
  
  tibble::tibble(
    Tumor_type = tumor,
    amplified_in_tumor = amplified_in_tumor,
    total_in_tumor = total_in_tumor,
    amplified_in_other = amplified_in_other,
    total_in_other = total_in_other,
    percent_in_tumor = 100 * amplified_in_tumor / total_in_tumor,
    percent_in_other = 100 * amplified_in_other / total_in_other,
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  )
}

# 3. LOAD DATA
rms_survival_all <- read.csv(
  survival_file,
  stringsAsFactors = FALSE
)

rms_survival_all <- rms_survival_all[, 1:8]

genes_altered <- read.delim(
  gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_gene_cols <- c(
  "studyID:sampleId",
  "Altered",
  gene_cols
)

missing_gene_cols <- setdiff(
  required_gene_cols,
  colnames(genes_altered)
)

if (length(missing_gene_cols) > 0) {
  stop(
    "Missing columns in gene alteration file: ",
    paste(missing_gene_cols, collapse = ", ")
  )
}

# 4. PATIENT-LEVEL ALTERATION DATA
alter_patient <- genes_altered |>
  dplyr::mutate(
    Patient.ID = stringr::str_extract(
      `studyID:sampleId`,
      "P-\\d+"
    )
  )

tumor_annot <- rms_survival_all |>
  dplyr::mutate(
    Patient.ID = stringr::str_extract(
      Patient.ID,
      "P-\\d+"
    )
  ) |>
  dplyr::filter(
    !is.na(Patient.ID),
    !is.na(Tumor_type),
    Tumor_type != ""
  ) |>
  dplyr::distinct(
    Patient.ID,
    Tumor_type
  )

# Check that each patient has only one tumor subtype.
tumor_annot_check <- tumor_annot |>
  dplyr::count(
    Patient.ID,
    name = "n_tumor_types"
  ) |>
  dplyr::filter(
    n_tumor_types > 1
  )

if (nrow(tumor_annot_check) > 0) {
  
  stop(
    "Some patients have more than one tumor subtype annotation. ",
    "See patients_with_multiple_tumor_types.csv."
  )
}

alter_annotated <- alter_patient |>
  dplyr::left_join(
    tumor_annot,
    by = "Patient.ID"
  )

patient_level <- alter_annotated |>
  dplyr::filter(
    !is.na(Patient.ID),
    !is.na(Tumor_type)
  ) |>
  dplyr::group_by(
    Patient.ID,
    Tumor_type
  ) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(gene_cols),
      ~ as.integer(any(.x == 1, na.rm = TRUE))
    ),
    Altered = as.integer(
      any(Altered == 1, na.rm = TRUE)
    ),
    n_samples = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Tumor_type = factor(
      Tumor_type,
      levels = tumor_order
    )
  )


patient_counts <- patient_level |>
  dplyr::count(
    Tumor_type,
    name = "Number_of_patients"
  )

# 5. RECURRENT GENE FREQUENCY BY SUBTYPE
freq_table <- patient_level |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(gene_cols),
    names_to = "Gene",
    values_to = "Altered_status"
  ) |>
  dplyr::group_by(
    Tumor_type,
    Gene
  ) |>
  dplyr::summarise(
    altered_patients = sum(
      Altered_status == 1,
      na.rm = TRUE
    ),
    total_patients = dplyr::n(),
    fraction_altered = altered_patients / total_patients,
    percent_altered = 100 * fraction_altered,
    .groups = "drop"
  )


# 6. GENE ENRICHMENT BY SUBTYPE
# One-vs-rest Fisher exact test + BH correction
tumor_types <- levels(
  droplevels(
    patient_level$Tumor_type
  )
)

enrichment_results <- tidyr::expand_grid(
  Gene = gene_cols,
  Tumor_type = tumor_types
) |>
  dplyr::mutate(
    result = purrr::map2(
      Gene,
      Tumor_type,
      ~ run_fisher_one_vs_rest(
        patient_level,
        .x,
        .y
      )
    )
  ) |>
  dplyr::select(
    result
  ) |>
  tidyr::unnest(
    result
  ) |>
  dplyr::mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    Direction = dplyr::case_when(
      p_adj_BH < 0.05 & odds_ratio > 1 ~ "Enriched",
      p_adj_BH < 0.05 & odds_ratio < 1 ~ "Depleted",
      TRUE ~ "Not significant"
    )
  ) |>
  dplyr::arrange(
    p_adj_BH,
    p_value
  )


# 7. GENE FREQUENCY HEATMAP
heatmap_gene_order <- c(
  "MYCN",
  "MDM2",
  "GLI1",
  "MYC",
  "CDK4",
  "TP53"
)

star_df <- enrichment_results |>
  dplyr::filter(
    Gene %in% heatmap_gene_order
  ) |>
  dplyr::mutate(
    sig_star = dplyr::case_when(
      odds_ratio > 1 & p_adj_BH < 0.001 ~ "***",
      odds_ratio > 1 & p_adj_BH < 0.01 ~ "**",
      odds_ratio > 1 & p_adj_BH < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  dplyr::select(
    Gene,
    Tumor_type,
    sig_star
  )

plot_df <- freq_table |>
  dplyr::filter(
    Gene %in% heatmap_gene_order
  ) |>
  dplyr::left_join(
    star_df,
    by = c(
      "Gene",
      "Tumor_type"
    )
  ) |>
  dplyr::mutate(
    Tumor_type = factor(
      Tumor_type,
      levels = tumor_order
    ),
    Gene = factor(
      Gene,
      levels = rev(
        heatmap_gene_order
      )
    ),
    label = dplyr::if_else(
      sig_star != "",
      paste0(
        sprintf(
          "%.1f%%",
          percent_altered
        ),
        "\n",
        sig_star
      ),
      sprintf(
        "%.1f%%",
        percent_altered
      )
    )
  )

heatmap_freq <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = Tumor_type,
    y = Gene,
    fill = percent_altered
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = label
    ),
    size = 4,
    lineheight = 0.85
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "red3",
    name = "% altered"
  ) +
  ggplot2::labs(
    x = "Tumor subtype",
    y = "Gene"
  ) +
  ggplot2::theme_classic(
    base_size = 12
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(
      color = "black",
      face = "italic"
    ),
    axis.title = ggplot2::element_text(
      color = "black"
    )
  )


# 8. 12q AMPLIFICATION
# Defined as alteration of MDM2, CDK4, or GLI1
patient_level_12q <- patient_level |>
  dplyr::mutate(
    Amp_12q = as.integer(
      dplyr::if_any(
        dplyr::all_of(genes_12q),
        ~ .x == 1
      )
    )
  )

freq_12q <- patient_level_12q |>
  dplyr::group_by(
    Tumor_type
  ) |>
  dplyr::summarise(
    amplified_patients = sum(
      Amp_12q == 1,
      na.rm = TRUE
    ),
    total_patients = dplyr::n(),
    fraction_amplified = amplified_patients / total_patients,
    percent_amplified = 100 * fraction_amplified,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Tumor_type = factor(
      Tumor_type,
      levels = tumor_order
    ),
    label = paste0(
      amplified_patients,
      "/",
      total_patients,
      "\n",
      sprintf(
        "%.1f%%",
        percent_amplified
      )
    )
  )


p_12q <- ggplot2::ggplot(
  freq_12q,
  ggplot2::aes(
    x = Tumor_type,
    y = fraction_amplified,
    fill = Tumor_type
  )
) +
  ggplot2::geom_col(
    width = 0.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = label
    ),
    vjust = -0.35,
    size = 4
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.15
      )
    )
  ) +
  ggplot2::scale_fill_manual(
    values = tumor_colors
  ) +
  ggplot2::labs(
    x = NULL,
    y = "12q amplified (%)"
  ) +
  ggplot2::theme_classic(
    base_size = 14
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(
      color = "black"
    ),
    axis.title.y = ggplot2::element_text(
      color = "black"
    )
  )


# 9. 12q ENRICHMENT BY SUBTYPE
enrichment_12q <- purrr::map_dfr(
  tumor_types,
  ~ run_fisher_12q(
    patient_level_12q,
    .x
  )
) |>
  dplyr::mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    Direction = dplyr::case_when(
      odds_ratio > 1 & p_adj_BH < 0.05 ~ "Enriched",
      odds_ratio < 1 & p_adj_BH < 0.05 ~ "Depleted",
      TRUE ~ "Not significant"
    )
  ) |>
  dplyr::arrange(
    p_adj_BH,
    p_value
  )


# 10. OVERALL ALTERATION FREQUENCY
overall_freq <- patient_level_12q |>
  dplyr::select(
    dplyr::all_of(gene_cols),
    Amp_12q
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Alteration",
    values_to = "Altered_status"
  ) |>
  dplyr::group_by(
    Alteration
  ) |>
  dplyr::summarise(
    altered_patients = sum(
      Altered_status == 1,
      na.rm = TRUE
    ),
    total_patients = dplyr::n(),
    percent_altered = 100 * altered_patients / total_patients,
    fraction = paste0(
      altered_patients,
      "/",
      total_patients
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Alteration = dplyr::recode(
      Alteration,
      Amp_12q = "12q amplification"
    )
  ) |>
  dplyr::arrange(
    dplyr::desc(
      percent_altered
    )
  )



# 12. SUBTYPE-SPECIFIC ONCOPRINTS
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(grid)
})


rms_survival_unique <- rms_survival_all |>
  dplyr::distinct(Patient.ID, .keep_all = TRUE)

# ARMS ONCOPRINT
ARMS_oncoprint <- read.csv(
  arms_oncoprint_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ARMS_oncoprint[ARMS_oncoprint == "no alteration"] <- ""
ARMS_oncoprint[ARMS_oncoprint == "not profiled"] <- ""

# Match your file's exact patient column name.
if ("Patient ID" %in% colnames(ARMS_oncoprint)) {
  colnames(ARMS_oncoprint)[
    colnames(ARMS_oncoprint) == "Patient ID"
  ] <- "Patient.ID"
}

ARMS_onco <- dplyr::left_join(
  ARMS_oncoprint,
  rms_survival_unique,
  by = "Patient.ID"
)

mat <- as.matrix(ARMS_onco)
mat <- t(ARMS_onco[4:13])

colnames(mat) <- ARMS_onco[, 3]

alter_fun <- list(
  background = function(x, y, w, h)
    grid::grid.rect(
      x, y, w, h,
      gp = grid::gpar(
        fill = "#f0f0f0",
        col = "white",
        lwd = 1.5
      )
    ),
  "AMP (driver)" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.85,
      gp = grid::gpar(fill = "#b2182b", col = NA)
    ),
  "missense" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "darkgreen", col = NA)
    ),
  "HOMDEL (driver)" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#436EEE", col = NA)
    ),
  "splice site" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#FF4500", col = NA)
    ),
  "nonsense" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.25,
      gp = grid::gpar(fill = "black", col = NA)
    )
)

col <- c(
  "AMP (driver)" = "#b2182b",
  "nonsense" = "black",
  "missense" = "darkgreen",
  "HOMDEL (driver)" = "#436EEE",
  "splice site" = "#FF4500"
)

partner <- ARMS_onco$PAX3.PAX7
partner[is.na(partner) | partner == ""] <- "Unknown"

status <- ARMS_onco$Status..NED..AWD.DOD.
status[is.na(status) | status == ""] <- "Unknown"

ha_top <- ComplexHeatmap::HeatmapAnnotation(
  Partner = partner,
  Status = status,
  col = list(
    Partner = c(
      "PAX3" = "#1f78b4",
      "PAX7" = "#33a02c",
      "Unknown" = "grey80"
    ),
    Status = c(
      "NED" = "#4daf4a",
      "AWD" = "#ff7f00",
      "DOD" = "#e41a1c",
      "Unknown" = "grey80"
    )
  )
)


ComplexHeatmap::oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  row_names_gp = grid::gpar(
    fontsize = 12,
    fontface = "italic"
  ),
  top_annotation = ha_top,
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE
)



# ERMS ONCOPRINT
erms_oncoprint <- read.csv(
  erms_oncoprint_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

erms_oncoprint[erms_oncoprint == "no alteration"] <- ""
erms_oncoprint[erms_oncoprint == "not profiled"] <- ""

if ("Patient ID" %in% colnames(erms_oncoprint)) {
  colnames(erms_oncoprint)[
    colnames(erms_oncoprint) == "Patient ID"
  ] <- "Patient.ID"
}

erms_oncoprint <- dplyr::left_join(
  erms_oncoprint,
  rms_survival_unique,
  by = "Patient.ID"
)

mat <- as.matrix(erms_oncoprint)
mat <- t(erms_oncoprint[4:24])

colnames(mat) <- erms_oncoprint[, 3]

alter_fun <- list(
  background = function(x, y, w, h)
    grid::grid.rect(
      x, y, w, h,
      gp = grid::gpar(
        fill = "#f0f0f0",
        col = "white",
        lwd = 1.5
      )
    ),
  "AMP (driver)" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.85,
      gp = grid::gpar(fill = "#b2182b", col = NA)
    ),
  "AMP" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.85,
      gp = grid::gpar(fill = "#b2182b", col = NA)
    ),
  "missense" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "darkgreen", col = NA)
    ),
  "HOMDEL (driver)" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#436EEE", col = NA)
    ),
  "HOMDEL" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#436EEE", col = NA)
    ),
  "splice site" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#FF4500", col = NA)
    ),
  "frameshift insertion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#7CCD7C", col = NA)
    ),
  "frameshift deletion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "blue3", col = NA)
    ),
  "Intragenic deletion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#ae017e", col = NA)
    ),
  "Inframe deletion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#CDCD00", col = NA)
    ),
  "Fusion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#9F79EE", col = NA)
    ),
  "nonsense" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "black", col = NA)
    ),
  "germline" = function(x, y, w, h) {
    grid::grid.segments(
      x0 = x - w * 0.42,
      y0 = y - h * 0.42,
      x1 = x + w * 0.42,
      y1 = y + h * 0.42,
      gp = grid::gpar(col = "black", lwd = 2)
    )
  }
)

col <- c(
  "AMP (driver)" = "#b2182b",
  "AMP" = "#b2182b",
  "missense" = "darkgreen",
  "nonsense" = "black",
  "HOMDEL (driver)" = "#436EEE",
  "HOMDEL" = "#436EEE",
  "splice site" = "#FF4500",
  "frameshift insertion" = "#7CCD7C",
  "frameshift deletion" = "blue3",
  "Inframe deletion" = "#CDCD00",
  "Intragenic deletion" = "#ae017e",
  "Fusion" = "#9F79EE",
  "germline" = "black"
)

status <- erms_oncoprint$Status..NED..AWD.DOD.
status[is.na(status) | status == ""] <- "Unknown"

ha_top <- ComplexHeatmap::HeatmapAnnotation(
  Status = status,
  col = list(
    Status = c(
      "NED" = "#4daf4a",
      "AWD" = "#ff7f00",
      "DOD" = "#e41a1c",
      "Unknown" = "grey80"
    )
  )
)

genes_order <- c(
  "TP53", "DICER1", "BCOR", "NF1", "FGFR4", "NRAS", "HRAS",
  "BRAF", "MDM2", "CDK4", "GLI1", "MYC", "MYCN", "CTNNB1",
  "FGFR1", "RAD21", "APC", "NSD3", "CDKN2A", "CDKN2B", "FBXW7"
)


ComplexHeatmap::oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  row_names_gp = grid::gpar(
    fontsize = 12,
    fontface = "italic"
  ),
  top_annotation = ha_top,
  row_order = genes_order,
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE
)



# SSRMS ONCOPRINT
ssrms_oncoprint <- read.csv(
  ssrms_oncoprint_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ssrms_oncoprint[ssrms_oncoprint == "no alteration"] <- ""
ssrms_oncoprint[ssrms_oncoprint == "not profiled"] <- ""

if ("Patient ID" %in% colnames(ssrms_oncoprint)) {
  colnames(ssrms_oncoprint)[
    colnames(ssrms_oncoprint) == "Patient ID"
  ] <- "Patient.ID"
}

ssrms_oncoprint <- dplyr::left_join(
  ssrms_oncoprint,
  rms_survival_unique,
  by = "Patient.ID"
)

mat <- as.matrix(ssrms_oncoprint)
mat <- t(ssrms_oncoprint[4:16])

colnames(mat) <- ssrms_oncoprint[, 3]

alter_fun <- list(
  background = function(x, y, w, h)
    grid::grid.rect(
      x, y, w, h,
      gp = grid::gpar(
        fill = "#f0f0f0",
        col = "white",
        lwd = 1.5
      )
    ),
  "AMP (driver)" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.85,
      gp = grid::gpar(fill = "#b2182b", col = NA)
    ),
  "AMP" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.85,
      gp = grid::gpar(fill = "#b2182b", col = NA)
    ),
  "missense" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "darkgreen", col = NA)
    ),
  "HOMDEL (driver)" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#436EEE", col = NA)
    ),
  "HOMDEL" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#436EEE", col = NA)
    ),
  "splice site" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.85, h * 0.55,
      gp = grid::gpar(fill = "#FF4500", col = NA)
    ),
  "frameshift insertion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#7CCD7C", col = NA)
    ),
  "L122R (driver)" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#7CC", col = NA)
    ),
  "frameshift deletion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "blue3", col = NA)
    ),
  "Intragenic deletion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#ae017e", col = NA)
    ),
  "Inframe deletion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#CDCD00", col = NA)
    ),
  "Fusion" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "#9F79EE", col = NA)
    ),
  "nonsense" = function(x, y, w, h)
    grid::grid.rect(
      x, y, w * 0.65, h * 0.25,
      gp = grid::gpar(fill = "black", col = NA)
    ),
  "germline" = function(x, y, w, h) {
    grid::grid.segments(
      x0 = x - w * 0.42,
      y0 = y - h * 0.42,
      x1 = x + w * 0.42,
      y1 = y + h * 0.42,
      gp = grid::gpar(col = "black", lwd = 2)
    )
  }
)

col <- c(
  "AMP (driver)" = "#b2182b",
  "AMP" = "#b2182b",
  "missense" = "darkgreen",
  "nonsense" = "black",
  "HOMDEL (driver)" = "#436EEE",
  "HOMDEL" = "#436EEE",
  "splice site" = "#FF4500",
  "frameshift insertion" = "#7CCD7C",
  "L122R (driver)" = "#7CC",
  "frameshift deletion" = "blue3",
  "Inframe deletion" = "#CDCD00",
  "Intragenic deletion" = "#ae017e",
  "Fusion" = "#9F79EE",
  "germline" = "black"
)

status <- ssrms_oncoprint$Status..NED..AWD.DOD.
status[is.na(status) | status == ""] <- "Unknown"

Driver <- ssrms_oncoprint$PAX3.PAX7
Driver[is.na(Driver) | Driver == ""] <- "Unknown"

driver_cols <- c(
  "MYOD1" = "#984ea3",
  "MEIS1-NCOA2" = "#377eb8",
  "FUS-TFCP2" = "#e41a1c",
  "ZFP64-NCOA3" = "#ff7f00",
  "Unknown" = "grey80"
)

ha_top <- ComplexHeatmap::HeatmapAnnotation(
  Status = status,
  Driver = Driver,
  col = list(
    Status = c(
      "NED" = "#4daf4a",
      "AWD" = "#ff7f00",
      "DOD" = "#e41a1c",
      "Unknown" = "grey80"
    ),
    Driver = driver_cols
  )
)

genes <- c(
  "MYOD1", "MDM2", "CDK4", "GLI1", "TP53",
  "CDKN2A", "CDKN2B", "MTAP", "ALK",
  "NRAS", "FGFR4", "MYC", "MYCN"
)


ComplexHeatmap::oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = col,
  row_names_gp = grid::gpar(
    fontsize = 12,
    fontface = "italic"
  ),
  top_annotation = ha_top,
  row_order = genes,
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE
)




# =============================================================================
# 13. SURVIVAL ANALYSIS — THREE AMPLIFICATION CLASSES
# No focal amplification vs chromosomal amplification vs ecDNA
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(survival)
  library(survminer)
  library(coxphf)
})


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

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA else x[1]
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

extract_cox_table <- function(fit, model_name) {
  s <- summary(fit)
  
  tibble::tibble(
    Model = model_name,
    Term = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"],
    CI_lower = s$conf.int[, "lower .95"],
    CI_upper = s$conf.int[, "upper .95"],
    P_value = s$coefficients[, "Pr(>|z|)"]
  )
}


# =============================================================================
# 13. SURVIVAL ANALYSIS — THREE AMPLIFICATION CLASSES
# No focal amplification vs chromosomal amplification vs ecDNA
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(survival)
  library(survminer)
  library(coxphf)
})


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

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA else x[1]
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

extract_cox_table <- function(fit, model_name) {
  s <- summary(fit)
  
  tibble::tibble(
    Model = model_name,
    Term = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"],
    CI_lower = s$conf.int[, "lower .95"],
    CI_upper = s$conf.int[, "upper .95"],
    P_value = s$coefficients[, "Pr(>|z|)"]
  )
}

# -----------------------------------------------------------------------------
# 13A. Prepare GCAP data and collapse to one row per patient
# Hierarchy: ecDNA/circular > chromosomal/noncircular > no focal
# -----------------------------------------------------------------------------

path_purity_map <- path_purity |>
  dplyr::transmute(
    sample = as.character(Sample.ID),
    path_purity = suppressWarnings(as.numeric(TumorPurity))
  ) |>
  dplyr::group_by(sample) |>
  dplyr::summarise(
    path_purity = first_nonmissing(path_purity),
    .groups = "drop"
  )

gcap_filtered <- gcap_results |>
  dplyr::mutate(
    sample = as.character(sample),
    class = tolower(trimws(as.character(class))),
    purity = suppressWarnings(as.numeric(purity))
  ) |>
  dplyr::left_join(
    path_purity_map,
    by = "sample"
  ) |>
  dplyr::filter(
    class == "circular" |
      (
        !is.na(purity) & purity >= 0.30 &
          !is.na(path_purity) & path_purity >= 30
      )
  )

clinical_gcap <- rms_survival_all |>
  dplyr::select(
    Sample.ID,
    Patient.ID,
    Tumor_type,
    Status..NED..AWD.DOD.,
    Time_of_FU
  ) |>
  dplyr::distinct()

gcap <- gcap_filtered |>
  dplyr::left_join(
    clinical_gcap,
    by = c("sample" = "Sample.ID")
  )


#remove samples with no clinical
gcap <- gcap |>
  dplyr::filter(!is.na(Patient.ID))

gcap_patient <- gcap |>
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
    Tumor_type = first_nonmissing(Tumor_type),
    Status = first_nonmissing(Status..NED..AWD.DOD.),
    Time_of_FU = safe_max(Time_of_FU),
    n_profiled_samples = dplyr::n(),
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
    GCAP_class = factor(
      class,
      levels = c("nofocal", "noncircular", "circular"),
      labels = c(
        "No focal amplification",
        "Chromosomal amplification",
        "ecDNA"
      )
    ),
    Tumor_type = factor(
      Tumor_type,
      levels = c("ARMS", "ERMS", "SSRMS")
    )
  )


print(table(gcap_patient$GCAP_class, useNA = "ifany"))
print(table(gcap_patient$Tumor_type, gcap_patient$GCAP_class, useNA = "ifany"))
print(table(gcap_patient$GCAP_class, gcap_patient$event))
print(table(gcap$class, useNA = "ifany"))

km_colors <- c(
  "No focal amplification" = "#008B45FF",
  "Chromosomal amplification" = "#EE0000FF",
  "ecDNA" = "#3B4992FF"
)

# -----------------------------------------------------------------------------
# 13B. Overall KM + global and raw pairwise log-rank tests
# -----------------------------------------------------------------------------

fit_3group_all <- survival::survfit(
  survival::Surv(Time_of_FU, event) ~ GCAP_class,
  data = gcap_patient
)

logrank_3group_all <- survival::survdiff(
  survival::Surv(Time_of_FU, event) ~ GCAP_class,
  data = gcap_patient
)

global_3group_p <- stats::pchisq(
  logrank_3group_all$chisq,
  df = length(logrank_3group_all$n) - 1,
  lower.tail = FALSE
)

global_3group <- tibble::tibble(
  Analysis = "Overall",
  Test = "Global log-rank",
  Chi_square = unname(logrank_3group_all$chisq),
  df = length(logrank_3group_all$n) - 1,
  P_value = global_3group_p
)

class_pairs <- combn(
  levels(gcap_patient$GCAP_class),
  2,
  simplify = FALSE
)

pairwise_3group <- purrr::map_dfr(
  class_pairs,
  function(pair) {
    
    df <- gcap_patient |>
      dplyr::filter(GCAP_class %in% pair) |>
      droplevels()
    
    lr <- survival::survdiff(
      survival::Surv(Time_of_FU, event) ~ GCAP_class,
      data = df
    )
    
    p <- stats::pchisq(
      lr$chisq,
      df = 1,
      lower.tail = FALSE
    )
    
    tibble::tibble(
      Tumor_type = "Overall",
      Group1 = pair[1],
      Group2 = pair[2],
      Group1_n = sum(df$GCAP_class == pair[1]),
      Group2_n = sum(df$GCAP_class == pair[2]),
      Group1_events = sum(
        df$event[df$GCAP_class == pair[1]],
        na.rm = TRUE
      ),
      Group2_events = sum(
        df$event[df$GCAP_class == pair[2]],
        na.rm = TRUE
      ),
      Chi_square = unname(lr$chisq),
      P_value = p
    )
  }
)



print(global_3group)
print(pairwise_3group)

pairwise_3group_labels <- pairwise_3group |>
  dplyr::mutate(
    label = paste0(
      Group1, " vs ", Group2,
      ": P = ", format_p(P_value)
    )
  )

pairwise_3group_text <- paste(
  pairwise_3group_labels$label,
  collapse = "\n"
)

p_all <- survminer::ggsurvplot(
  fit_3group_all,
  data = gcap_patient,
  pval = paste0(
    "Global log-rank P = ",
    format_p(global_3group_p)
  ),
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = T,
  xlab = "Time since diagnosis (months)",
  ylab = "Overall survival probability",
  palette = unname(km_colors),
  legend.title = "Amplification class",
  legend.labs = names(km_colors),
  legend = "right",
  censor.shape = 124,
  censor.size = 3,
  break.time.by = 24,
  ggtheme = ggplot2::theme_classic(base_size = 13),
  tables.theme = ggplot2::theme_classic(base_size = 10),
  title = ""
)

p_all$plot <- p_all$plot +
  ggplot2::annotate(
    "text",
    x = max(
      gcap_patient$Time_of_FU,
      na.rm = TRUE
    ) * 0.48,
    y = 0.22,
    label = paste0(
      "Pairwise log-rank:\n",
      pairwise_3group_text
    ),
    hjust = 0,
    size = 3.0
  )

print(p_all)

print(p_all)

# -----------------------------------------------------------------------------
# 13C. Overall Cox models: unadjusted and adjusted for RMS subtype
# -----------------------------------------------------------------------------

cox_3group_unadjusted <- survival::coxph(
  survival::Surv(Time_of_FU, event) ~ GCAP_class,
  data = gcap_patient,
  ties = "efron"
)

cox_3group_adjusted_reduced <- survival::coxph(
  survival::Surv(Time_of_FU, event) ~ Tumor_type,
  data = gcap_patient,
  ties = "efron"
)

cox_3group_adjusted <- survival::coxph(
  survival::Surv(Time_of_FU, event) ~ GCAP_class + Tumor_type,
  data = gcap_patient,
  ties = "efron"
)

cox_3group_results <- dplyr::bind_rows(
  extract_cox_table(
    cox_3group_unadjusted,
    "Unadjusted"
  ),
  extract_cox_table(
    cox_3group_adjusted,
    "Adjusted for RMS subtype"
  )
)

cox_3group_global <- tibble::tibble(
  Model = c(
    "Unadjusted",
    "Adjusted for RMS subtype"
  ),
  Comparison = "GCAP class overall",
  P_value_LRT = c(
    stats::anova(
      survival::coxph(
        survival::Surv(Time_of_FU, event) ~ 1,
        data = gcap_patient
      ),
      cox_3group_unadjusted,
      test = "LRT"
    )[2, "Pr(>|Chi|)"],
    stats::anova(
      cox_3group_adjusted_reduced,
      cox_3group_adjusted,
      test = "LRT"
    )[2, "Pr(>|Chi|)"]
  )
)

ph_3group <- dplyr::bind_rows(
  as.data.frame(
    survival::cox.zph(cox_3group_unadjusted)$table
  ) |>
    tibble::rownames_to_column("Term") |>
    dplyr::mutate(Model = "Unadjusted", .before = 1),
  
  as.data.frame(
    survival::cox.zph(cox_3group_adjusted)$table
  ) |>
    tibble::rownames_to_column("Term") |>
    dplyr::mutate(Model = "Adjusted for RMS subtype", .before = 1)
)




print(cox_3group_results)
print(cox_3group_global)
print(ph_3group)

# -----------------------------------------------------------------------------
# 13D. Global and pairwise log-rank tests within each RMS subtype
# -----------------------------------------------------------------------------

tumor_types_surv <- levels(
  droplevels(gcap_patient$Tumor_type)
)

subtype_global_3group <- purrr::map_dfr(
  tumor_types_surv,
  function(tumor) {
    
    df <- gcap_patient |>
      dplyr::filter(Tumor_type == tumor) |>
      droplevels()
    
    if (dplyr::n_distinct(df$GCAP_class) < 2) {
      return(tibble::tibble())
    }
    
    lr <- survival::survdiff(
      survival::Surv(Time_of_FU, event) ~ GCAP_class,
      data = df
    )
    
    tibble::tibble(
      Tumor_type = tumor,
      n = nrow(df),
      events = sum(df$event, na.rm = TRUE),
      Chi_square = unname(lr$chisq),
      df = length(lr$n) - 1,
      P_value = stats::pchisq(
        lr$chisq,
        df = length(lr$n) - 1,
        lower.tail = FALSE
      )
    )
  }
)

subtype_pairwise_3group <- purrr::map_dfr(
  tumor_types_surv,
  function(tumor) {
    
    df <- gcap_patient |>
      dplyr::filter(Tumor_type == tumor) |>
      droplevels()
    
    classes_here <- levels(droplevels(df$GCAP_class))
    
    if (length(classes_here) < 2) {
      return(tibble::tibble())
    }
    
    pairs_here <- combn(
      classes_here,
      2,
      simplify = FALSE
    )
    
    purrr::map_dfr(
      pairs_here,
      function(pair) {
        
        df_pair <- df |>
          dplyr::filter(GCAP_class %in% pair) |>
          droplevels()
        
        lr <- survival::survdiff(
          survival::Surv(Time_of_FU, event) ~ GCAP_class,
          data = df_pair
        )
        
        tibble::tibble(
          Tumor_type = tumor,
          Group1 = pair[1],
          Group2 = pair[2],
          Group1_n = sum(df_pair$GCAP_class == pair[1]),
          Group2_n = sum(df_pair$GCAP_class == pair[2]),
          Group1_events = sum(
            df_pair$event[df_pair$GCAP_class == pair[1]],
            na.rm = TRUE
          ),
          Group2_events = sum(
            df_pair$event[df_pair$GCAP_class == pair[2]],
            na.rm = TRUE
          ),
          P_value = stats::pchisq(
            lr$chisq,
            df = 1,
            lower.tail = FALSE
          )
        )
      }
    )
  }
)


# FDR correction specifically within ERMS across the three pairwise
# three-group comparisons.
erms_pairwise_3group_fdr <- subtype_pairwise_3group |>
  dplyr::filter(Tumor_type == "ERMS") |>
  dplyr::mutate(
    FDR = stats::p.adjust(P_value, method = "BH"),
    significance = dplyr::case_when(
      FDR <= 0.001 ~ "***",
      FDR <= 0.01 ~ "**",
      FDR <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  dplyr::arrange(FDR, P_value)



print(subtype_global_3group)
print(subtype_pairwise_3group)
print(erms_pairwise_3group_fdr)


# Subtype-specific KM figures.
plots_3group_by_type <- lapply(
  tumor_types_surv,
  function(tumor) {
    
    df <- gcap_patient |>
      dplyr::filter(
        Tumor_type == tumor,
        !is.na(GCAP_class),
        !is.na(Time_of_FU),
        !is.na(event)
      ) |>
      droplevels()
    
    if (dplyr::n_distinct(df$GCAP_class) < 2) {
      return(NULL)
    }
    
    fit <- survival::survfit(
      survival::Surv(
        Time_of_FU,
        event
      ) ~ GCAP_class,
      data = df
    )
    
    lvls <- levels(df$GCAP_class)
    
    p <- survminer::ggsurvplot(
      fit,
      data = df,
      
      pval = TRUE,
      conf.int = FALSE,
      
      risk.table = TRUE,
      risk.table.height = 0.25,
      risk.table.y.text = TRUE,
      risk.table.y.text.col = TRUE,
      
      xlab = "Time since diagnosis (months)",
      ylab = "Overall survival probability",
      
      palette = unname(
        km_colors[lvls]
      ),
      
      legend.title = "Amplification class",
      legend.labs = lvls,
      legend = "right",
      
      censor = TRUE,
      censor.shape = 3,
      censor.size = 2.5,
      
      break.time.by = 24,
      
      ggtheme = ggplot2::theme_classic(
        base_size = 13
      ),
      
      tables.theme = ggplot2::theme_classic(
        base_size = 10
      ),
      
      title = tumor
    )
    
    
    # Match formatting of the overall KM.
    p$plot <- p$plot +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 14,
          hjust = 0.5
        ),
        legend.position = "right",
        legend.text = ggplot2::element_text(
          size = 10
        ),
        legend.title = ggplot2::element_text(
          size = 11
        ),
        axis.title = ggplot2::element_text(
          size = 12
        ),
        axis.text = ggplot2::element_text(
          size = 10
        )
      )
    
    
    # ERMS pairwise log-rank results.
    if (tumor == "ERMS") {
      
      pairwise_erms <- erms_pairwise_3group_fdr |>
        dplyr::mutate(
          label = paste0(
            Group1,
            " vs ",
            Group2,
            ": FDR = ",
            format_p(FDR),
            " ",
            significance
          )
        )
      
      if (nrow(pairwise_erms) > 0) {
        
        p$plot <- p$plot +
          ggplot2::annotate(
            "text",
            x = max(
              df$Time_of_FU,
              na.rm = TRUE
            ) * 0.48,
            y = 0.22,
            label = paste0(
              "Pairwise log-rank:\n",
              paste(
                pairwise_erms$label,
                collapse = "\n"
              )
            ),
            hjust = 0,
            size = 3.0
          )
      }
    }
    
    p
  }
)

names(plots_3group_by_type) <- tumor_types_surv

plots_3group_by_type <- plots_3group_by_type[
  !vapply(
    plots_3group_by_type,
    is.null,
    logical(1)
  )
]


# ============================================================
# Combined subtype figure
# ============================================================

if (length(plots_3group_by_type) > 0) {
  
  p_3group_split <- survminer::arrange_ggsurvplots(
    plots_3group_by_type,
    print = FALSE,
    ncol = 1,
    nrow = length(plots_3group_by_type)
  )
  
  
  print(p_3group_split)
  
  
  
  # Multi-page PDF: one subtype per page.
  
  for (tumor in names(plots_3group_by_type)) {
    print(
      plots_3group_by_type[[tumor]]
    )
  }
  
}
# -----------------------------------------------------------------------------
# 13E. Pairwise Firth Cox models overall and within each subtype
# -----------------------------------------------------------------------------

run_pairwise_firth_3group <- function(df, tumor_label, pair) {
  
  df_pair <- df |>
    dplyr::filter(GCAP_class %in% pair) |>
    droplevels()
  
  if (dplyr::n_distinct(df_pair$GCAP_class) < 2) {
    return(NULL)
  }
  
  df_pair$GCAP_class <- stats::relevel(
    df_pair$GCAP_class,
    ref = pair[1]
  )
  
  fit <- tryCatch(
    coxphf::coxphf(
      survival::Surv(Time_of_FU, event) ~ GCAP_class,
      data = df_pair
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(NULL)
  }
  
  tibble::tibble(
    Tumor_type = tumor_label,
    Reference = pair[1],
    Comparison = pair[2],
    n_reference = sum(df_pair$GCAP_class == pair[1]),
    n_comparison = sum(df_pair$GCAP_class == pair[2]),
    events_reference = sum(
      df_pair$event[df_pair$GCAP_class == pair[1]],
      na.rm = TRUE
    ),
    events_comparison = sum(
      df_pair$event[df_pair$GCAP_class == pair[2]],
      na.rm = TRUE
    ),
    HR_Firth = exp(fit$coefficients[1]),
    CI_lower = fit$ci.lower[1],
    CI_upper = fit$ci.upper[1],
    P_value = fit$prob[1]
  )
}

pairwise_firth_3group <- purrr::map_dfr(
  class_pairs,
  ~ run_pairwise_firth_3group(
    gcap_patient,
    "Overall",
    .x
  )
)

for (tumor in tumor_types_surv) {
  
  df <- gcap_patient |>
    dplyr::filter(Tumor_type == tumor) |>
    droplevels()
  
  classes_here <- levels(droplevels(df$GCAP_class))
  
  if (length(classes_here) < 2) {
    next
  }
  
  pairs_here <- combn(
    classes_here,
    2,
    simplify = FALSE
  )
  
  tmp <- purrr::map_dfr(
    pairs_here,
    ~ run_pairwise_firth_3group(
      df,
      tumor,
      .x
    )
  )
  
  pairwise_firth_3group <- dplyr::bind_rows(
    pairwise_firth_3group,
    tmp
  )
}

erms_pairwise_firth_fdr <- pairwise_firth_3group |>
  dplyr::filter(Tumor_type == "ERMS") |>
  dplyr::mutate(
    FDR = stats::p.adjust(P_value, method = "BH"),
    significance = dplyr::case_when(
      FDR <= 0.001 ~ "***",
      FDR <= 0.01 ~ "**",
      FDR <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  dplyr::arrange(FDR, P_value)



print(pairwise_firth_3group)
print(erms_pairwise_firth_fdr)

# -----------------------------------------------------------------------------

# =============================================================================
# 14. SURVIVAL ANALYSIS — BINARY ecDNA VS NON-ecDNA
# non-ecDNA combines no focal amplification + chromosomal amplification
# The same survival tests are repeated using this binary definition.
# =============================================================================

gcap_patient_binary <- gcap_patient |>
  dplyr::mutate(
    ecDNA_binary = factor(
      ifelse(class == "circular", "ecDNA", "non-ecDNA"),
      levels = c("non-ecDNA", "ecDNA")
    )
  )

binary_colors <- c(
  "non-ecDNA" = "#7F7F7F",
  "ecDNA" = "#3B4992FF"
)

# -----------------------------------------------------------------------------
# 14A. Overall binary KM + log-rank
# -----------------------------------------------------------------------------

fit_binary_all <- survival::survfit(
  survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
  data = gcap_patient_binary
)

logrank_binary_all <- survival::survdiff(
  survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
  data = gcap_patient_binary
)

binary_logrank_p <- stats::pchisq(
  logrank_binary_all$chisq,
  df = 1,
  lower.tail = FALSE
)

binary_overall_logrank <- tibble::tibble(
  Analysis = "Overall ecDNA vs non-ecDNA",
  n_non_ecDNA = sum(
    gcap_patient_binary$ecDNA_binary == "non-ecDNA"
  ),
  n_ecDNA = sum(
    gcap_patient_binary$ecDNA_binary == "ecDNA"
  ),
  events_non_ecDNA = sum(
    gcap_patient_binary$event[
      gcap_patient_binary$ecDNA_binary == "non-ecDNA"
    ],
    na.rm = TRUE
  ),
  events_ecDNA = sum(
    gcap_patient_binary$event[
      gcap_patient_binary$ecDNA_binary == "ecDNA"
    ],
    na.rm = TRUE
  ),
  Chi_square = unname(logrank_binary_all$chisq),
  P_value = binary_logrank_p
)


print(binary_overall_logrank)

p_binary <- survminer::ggsurvplot(
  fit_binary_all,
  data = gcap_patient_binary,
  pval = paste0(
    "Log-rank P = ",
    format_p(binary_logrank_p)
  ),
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  xlab = "Time since diagnosis (months)",
  ylab = "Overall survival probability",
  palette = unname(binary_colors),
  legend.title = "ecDNA status",
  legend.labs = names(binary_colors),
  legend = "right",
  censor.shape = 124,
  censor.size = 3,
  break.time.by = 24,
  ggtheme = ggplot2::theme_classic(base_size = 13),
  tables.theme = ggplot2::theme_classic(base_size = 10),
  title = ""
)

print(p_binary)

print(p_binary)

# -----------------------------------------------------------------------------
# 14B. Overall binary Cox models: unadjusted and adjusted for RMS subtype
# -----------------------------------------------------------------------------

cox_binary_unadjusted <- survival::coxph(
  survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
  data = gcap_patient_binary,
  ties = "efron"
)

cox_binary_adjusted <- survival::coxph(
  survival::Surv(Time_of_FU, event) ~ ecDNA_binary + Tumor_type,
  data = gcap_patient_binary,
  ties = "efron"
)

cox_binary_results <- dplyr::bind_rows(
  extract_cox_table(
    cox_binary_unadjusted,
    "Unadjusted"
  ),
  extract_cox_table(
    cox_binary_adjusted,
    "Adjusted for RMS subtype"
  )
)

ph_binary <- dplyr::bind_rows(
  as.data.frame(
    survival::cox.zph(cox_binary_unadjusted)$table
  ) |>
    tibble::rownames_to_column("Term") |>
    dplyr::mutate(Model = "Unadjusted", .before = 1),
  
  as.data.frame(
    survival::cox.zph(cox_binary_adjusted)$table
  ) |>
    tibble::rownames_to_column("Term") |>
    dplyr::mutate(Model = "Adjusted for RMS subtype", .before = 1)
)



print(cox_binary_results)
print(ph_binary)

# Firth binary Cox overall.
firth_binary_overall <- coxphf::coxphf(
  survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
  data = gcap_patient_binary
)

firth_binary_overall_results <- tibble::tibble(
  Analysis = "Overall",
  HR_Firth = exp(firth_binary_overall$coefficients[1]),
  CI_lower = firth_binary_overall$ci.lower[1],
  CI_upper = firth_binary_overall$ci.upper[1],
  P_value = firth_binary_overall$prob[1]
)


print(firth_binary_overall_results)

# -----------------------------------------------------------------------------
# 14C. Binary ecDNA vs non-ecDNA within each RMS subtype
# Repeat log-rank + Firth Cox
# -----------------------------------------------------------------------------

binary_by_subtype <- purrr::map_dfr(
  tumor_types_surv,
  function(tumor) {
    
    df <- gcap_patient_binary |>
      dplyr::filter(Tumor_type == tumor) |>
      droplevels()
    
    if (dplyr::n_distinct(df$ecDNA_binary) < 2) {
      return(tibble::tibble())
    }
    
    lr <- survival::survdiff(
      survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
      data = df
    )
    
    p_lr <- stats::pchisq(
      lr$chisq,
      df = 1,
      lower.tail = FALSE
    )
    
    firth <- tryCatch(
      coxphf::coxphf(
        survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
        data = df
      ),
      error = function(e) NULL
    )
    
    tibble::tibble(
      Tumor_type = tumor,
      n_non_ecDNA = sum(df$ecDNA_binary == "non-ecDNA"),
      n_ecDNA = sum(df$ecDNA_binary == "ecDNA"),
      events_non_ecDNA = sum(
        df$event[df$ecDNA_binary == "non-ecDNA"],
        na.rm = TRUE
      ),
      events_ecDNA = sum(
        df$event[df$ecDNA_binary == "ecDNA"],
        na.rm = TRUE
      ),
      Logrank_P = p_lr,
      HR_Firth = if (is.null(firth)) NA_real_ else exp(firth$coefficients[1]),
      CI_lower = if (is.null(firth)) NA_real_ else firth$ci.lower[1],
      CI_upper = if (is.null(firth)) NA_real_ else firth$ci.upper[1],
      Firth_P = if (is.null(firth)) NA_real_ else firth$prob[1]
    )
  }
)


print(binary_by_subtype)

binary_plots_by_type <- lapply(
  tumor_types_surv,
  function(tumor) {
    
    df <- gcap_patient_binary |>
      dplyr::filter(Tumor_type == tumor) |>
      droplevels()
    
    if (dplyr::n_distinct(df$ecDNA_binary) < 2) {
      return(NULL)
    }
    
    fit <- survival::survfit(
      survival::Surv(Time_of_FU, event) ~ ecDNA_binary,
      data = df
    )
    
    survminer::ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      conf.int = FALSE,
      risk.table = FALSE,
      xlab = "Time since diagnosis (months)",
      ylab = "Overall survival probability",
      palette = unname(binary_colors),
      legend.title = "ecDNA status",
      legend.labs = names(binary_colors),
      censor.shape = 124,
      censor.size = 3,
      legend = "right",
      break.time.by = 24,
      ggtheme = ggplot2::theme_classic(base_size = 13),
      title = paste("Tumor type:", tumor)
    )
  }
)

binary_plots_by_type <- binary_plots_by_type[
  !vapply(binary_plots_by_type, is.null, logical(1))
]

if (length(binary_plots_by_type) > 0) {
  
  p_binary_split <- survminer::arrange_ggsurvplots(
    binary_plots_by_type,
    print = FALSE,
    ncol = 1,
    nrow = length(binary_plots_by_type)
  )
  
  print(p_binary_split)
}

# -----------------------------------------------------------------------------

# =============================================================================
# 15. ERMS MULTIVARIABLE CLINICAL MODEL + PUBLICATION-STYLE FOREST PLOT
#
# THREE amplification classes:
#   No focal amplification (reference)
#   Chromosomal amplification
#   ecDNA
#
# Covariates:
#   TP53 status
#   Age at diagnosis (<10 years vs >=10 years)
#   Sex
#   Local stage (Stage 1-2 vs Stage 3-4)
#   Clinical group (Group I-II vs Group III-IV)
#
# Germline status is NOT included in this model.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(survival)
  library(coxphf)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})


erms_data <- read.csv(
  erms_clinical_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(erms_data) <- make.names(
  names(erms_data),
  unique = TRUE
)

required_erms_clinical_cols <- c(
  "Patient.ID",
  "Local.Stage.Clinical.Group",
  "age.sex"
)

missing_erms_clinical_cols <- setdiff(
  required_erms_clinical_cols,
  names(erms_data)
)

if (length(missing_erms_clinical_cols) > 0) {
  stop(
    "Missing required ERMS clinical columns: ",
    paste(missing_erms_clinical_cols, collapse = ", "),
    "\nAvailable columns are: ",
    paste(names(erms_data), collapse = ", ")
  )
}

erms_data_clean <- erms_data[
  ,
  required_erms_clinical_cols,
  drop = FALSE
]

# -----------------------------------------------------------------------------
# 15A. Extract local stage, clinical group, age, and sex
# -----------------------------------------------------------------------------

stage_group_split <- stringr::str_split_fixed(
  trimws(
    as.character(
      erms_data_clean$Local.Stage.Clinical.Group
    )
  ),
  "/",
  2
)

age_sex_split <- stringr::str_split_fixed(
  trimws(
    as.character(
      erms_data_clean$age.sex
    )
  ),
  "/",
  2
)

erms_clinical_raw <- erms_data_clean |>
  dplyr::transmute(
    Patient.ID = trimws(
      as.character(Patient.ID)
    ),
    
    Local_stage = suppressWarnings(
      as.numeric(
        trimws(
          stage_group_split[, 1]
        )
      )
    ),
    
    Clinical_group = trimws(
      stage_group_split[, 2]
    ),
    
    Age = suppressWarnings(
      as.numeric(
        trimws(
          age_sex_split[, 1]
        )
      )
    ),
    
    Sex_code = toupper(
      trimws(
        age_sex_split[, 2]
      )
    )
  ) |>
  dplyr::mutate(
    Clinical_group = dplyr::na_if(
      Clinical_group,
      ""
    ),
    
    Sex_code = dplyr::na_if(
      Sex_code,
      ""
    )
  )

first_nonmissing_numeric <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    NA_real_
  } else {
    x[1]
  }
}

first_nonmissing_character <- function(x) {
  x <- as.character(x)
  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]
  
  if (length(x) == 0) {
    NA_character_
  } else {
    x[1]
  }
}

erms_clinical_covariates <- erms_clinical_raw |>
  dplyr::filter(
    !is.na(Patient.ID),
    Patient.ID != ""
  ) |>
  dplyr::group_by(
    Patient.ID
  ) |>
  dplyr::summarise(
    Local_stage = first_nonmissing_numeric(
      Local_stage
    ),
    
    Clinical_group = first_nonmissing_character(
      Clinical_group
    ),
    
    Age = first_nonmissing_numeric(
      Age
    ),
    
    Sex_code = first_nonmissing_character(
      Sex_code
    ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Clinical_group = factor(
      Clinical_group,
      levels = c(
        "I",
        "II",
        "III",
        "IV"
      )
    ),
    
    Clinical_group_num = dplyr::case_when(
      Clinical_group == "I" ~ 1,
      Clinical_group == "II" ~ 2,
      Clinical_group == "III" ~ 3,
      Clinical_group == "IV" ~ 4,
      TRUE ~ NA_real_
    ),
    
    Sex = factor(
      Sex_code,
      levels = c(
        "F",
        "M"
      ),
      labels = c(
        "Female",
        "Male"
      )
    ),
    
    # Age modeled categorically:
    #   <10 years = reference
    #   >=10 years = comparison
    Age_group = factor(
      dplyr::case_when(
        is.na(Age) ~ NA_character_,
        Age < 10 ~ "<10 years",
        Age >= 10 ~ ">=10 years"
      ),
      levels = c(
        "<10 years",
        ">=10 years"
      )
    )
  )


# -----------------------------------------------------------------------------
# 15B. Join clinical covariates + TP53 to THREE-GROUP GCAP data
# -----------------------------------------------------------------------------

erms_multivariable_3group_all <- gcap_patient |>
  dplyr::filter(
    Tumor_type == "ERMS"
  ) |>
  dplyr::left_join(
    patient_level |>
      dplyr::select(
        Patient.ID,
        TP53
      ),
    by = "Patient.ID"
  ) |>
  dplyr::left_join(
    erms_clinical_covariates,
    by = "Patient.ID"
  ) |>
  dplyr::mutate(
    TP53_status = factor(
      TP53,
      levels = c(
        0,
        1
      ),
      labels = c(
        "TP53 WT",
        "TP53 altered"
      )
    ),
    
    GCAP_class = factor(
      GCAP_class,
      levels = c(
        "No focal amplification",
        "Chromosomal amplification",
        "ecDNA"
      )
    ),
    
    Sex = factor(
      Sex,
      levels = c(
        "Female",
        "Male"
      )
    ),
    
    # Collapsed local stage:
    #   Stage 1-2 = reference
    #   Stage 3-4 = higher-stage group
    Stage_group = factor(
      dplyr::case_when(
        Local_stage %in% c(1, 2) ~ "Stage 1-2",
        Local_stage %in% c(3, 4) ~ "Stage 3-4",
        TRUE ~ NA_character_
      ),
      levels = c(
        "Stage 1-2",
        "Stage 3-4"
      )
    ),
    
    # Collapsed clinical group:
    #   Group I-II = reference
    #   Group III-IV = higher clinical-group category
    Clinical_group_binary = factor(
      dplyr::case_when(
        Clinical_group %in% c("I", "II") ~ "Group I-II",
        Clinical_group %in% c("III", "IV") ~ "Group III-IV",
        TRUE ~ NA_character_
      ),
      levels = c(
        "Group I-II",
        "Group III-IV"
      )
    )
  )

cat(
  "\nAll ERMS patients before clinical complete-case filtering: n=",
  nrow(
    erms_multivariable_3group_all
  ),
  "\n",
  sep = ""
)

# Missingness table.
erms_multivariable_missingness <- tibble::tibble(
  Variable = c(
    "GCAP class",
    "TP53",
    "Stage group",
    "Clinical group",
    "Age",
    "Sex",
    "Follow-up",
    "Event"
  ),
  
  Missing_n = c(
    sum(is.na(
      erms_multivariable_3group_all$GCAP_class
    )),
    sum(is.na(
      erms_multivariable_3group_all$TP53_status
    )),
    sum(is.na(
      erms_multivariable_3group_all$Stage_group
    )),
    sum(is.na(
      erms_multivariable_3group_all$Clinical_group_binary
    )),
    sum(is.na(
      erms_multivariable_3group_all$Age_group
    )),
    sum(is.na(
      erms_multivariable_3group_all$Sex
    )),
    sum(is.na(
      erms_multivariable_3group_all$Time_of_FU
    )),
    sum(is.na(
      erms_multivariable_3group_all$event
    ))
  )
) |>
  dplyr::mutate(
    Total_n = nrow(
      erms_multivariable_3group_all
    ),
    Available_n =
      Total_n -
      Missing_n
  )

print(
  erms_multivariable_missingness
)


# Patients missing any variable actually used in THIS model.
erms_missing_model_covariates <- erms_multivariable_3group_all |>
  dplyr::filter(
    is.na(GCAP_class) |
      is.na(TP53_status) |
      is.na(Stage_group) |
      is.na(Clinical_group_binary) |
      is.na(Age_group) |
      is.na(Sex) |
      is.na(Time_of_FU) |
      is.na(event)
  ) |>
  dplyr::select(
    Patient.ID,
    GCAP_class,
    TP53_status,
    Local_stage,
    Stage_group,
    Clinical_group,
    Clinical_group_binary,
    Age,
    Sex,
    Time_of_FU,
    event
  )


print(
  erms_missing_model_covariates
)

# Complete-case cohort ONLY for variables in the no-germline model.
erms_multivariable_3group <- erms_multivariable_3group_all |>
  dplyr::filter(
    !is.na(GCAP_class),
    !is.na(TP53_status),
    !is.na(Stage_group),
    !is.na(Clinical_group_binary),
    !is.na(Age_group),
    !is.na(Sex),
    !is.na(Time_of_FU),
    !is.na(event)
  ) |>
  droplevels()


cat(
  "\nERMS no-germline multivariable cohort: n=",
  nrow(
    erms_multivariable_3group
  ),
  "; deaths=",
  sum(
    erms_multivariable_3group$event,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "\nAge group by event in the model cohort:\n"
)

print(
  table(
    erms_multivariable_3group$Age_group,
    erms_multivariable_3group$event,
    useNA = "ifany"
  )
)


cat(
  "\nCollapsed stage by event:\n"
)

print(
  table(
    erms_multivariable_3group$Stage_group,
    erms_multivariable_3group$event,
    useNA = "ifany"
  )
)

cat(
  "\nCollapsed clinical group by event:\n"
)

print(
  table(
    erms_multivariable_3group$Clinical_group_binary,
    erms_multivariable_3group$event,
    useNA = "ifany"
  )
)

# -----------------------------------------------------------------------------
# 15C. PRIMARY Firth Cox model
# -----------------------------------------------------------------------------

firth_erms_3group_clinical <- coxphf::coxphf(
  survival::Surv(
    Time_of_FU,
    event
  ) ~
    GCAP_class +
    TP53_status +
    Age_group +
    Sex +
    Stage_group +
    Clinical_group_binary,
  
  data = erms_multivariable_3group,
  
  maxit = 1000,
  maxstep = 0.5
)

print(
  summary(
    firth_erms_3group_clinical
  )
)

erms_firth_3group_clinical_results <- tibble::tibble(
  Variable = names(
    firth_erms_3group_clinical$coefficients
  ),
  
  HR = exp(
    firth_erms_3group_clinical$coefficients
  ),
  
  CI_lower =
    firth_erms_3group_clinical$ci.lower,
  
  CI_upper =
    firth_erms_3group_clinical$ci.upper,
  
  P_value =
    firth_erms_3group_clinical$prob
)


print(
  erms_firth_3group_clinical_results
)

# -----------------------------------------------------------------------------
# 15D. Standard Cox for global GCAP LRT and PH diagnostics
# -----------------------------------------------------------------------------

cox_erms_3group_clinical_reduced <- survival::coxph(
  survival::Surv(
    Time_of_FU,
    event
  ) ~
    TP53_status +
    Age_group +
    Sex +
    Stage_group +
    Clinical_group_binary,
  
  data = erms_multivariable_3group,
  
  ties = "efron"
)

cox_erms_3group_clinical <- survival::coxph(
  survival::Surv(
    Time_of_FU,
    event
  ) ~
    GCAP_class +
    TP53_status +
    Age_group +
    Sex +
    Stage_group +
    Clinical_group_binary,
  
  data = erms_multivariable_3group,
  
  ties = "efron"
)

erms_3group_clinical_global_p <- stats::anova(
  cox_erms_3group_clinical_reduced,
  cox_erms_3group_clinical,
  test = "LRT"
)[2, "Pr(>|Chi|)"]

erms_standard_3group_clinical_results <- extract_cox_table(
  cox_erms_3group_clinical,
  "ERMS 3-group clinical multivariable Cox"
)

ph_erms_3group_clinical <- as.data.frame(
  survival::cox.zph(
    cox_erms_3group_clinical
  )$table
) |>
  tibble::rownames_to_column(
    "Term"
  )

erms_3group_global_adjusted <- tibble::tibble(
  Model =
    "ERMS clinical multivariable Cox",
  
  Comparison = paste0(
    "GCAP class overall adjusted for ",
    "TP53, age, sex, collapsed local stage, ",
    "and collapsed clinical group"
  ),
  
  P_value_LRT =
    erms_3group_clinical_global_p
)




print(
  erms_standard_3group_clinical_results
)

print(
  erms_3group_global_adjusted
)

print(
  ph_erms_3group_clinical
)

# =============================================================================
# 15E. PUBLICATION-STYLE FOREST PLOT
# =============================================================================

get_firth_term <- function(term_name) {
  
  idx <- match(
    term_name,
    erms_firth_3group_clinical_results$Variable
  )
  
  if (is.na(idx)) {
    return(
      tibble::tibble(
        HR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_
      )
    )
  }
  
  erms_firth_3group_clinical_results[
    idx,
    c(
      "HR",
      "CI_lower",
      "CI_upper",
      "P_value"
    )
  ]
}

fmt_hr_number <- function(x) {
  # Show up to 2 decimals, but remove unnecessary trailing zeros.
  # Examples: 5.50 -> 5.5; 1.30 -> 1.3; 1.35 -> 1.35
  sub(
    "\\.?0+$",
    "",
    sprintf("%.2f", x)
  )
}

fmt_hr_ci <- function(hr, lo, hi) {
  
  if (
    is.na(hr) ||
    is.na(lo) ||
    is.na(hi)
  ) {
    return("not estimable")
  }
  
  paste0(
    fmt_hr_number(hr),
    " (",
    fmt_hr_number(lo),
    "-",
    fmt_hr_number(hi),
    ")"
  )
}

fmt_model_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NE",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

sig_symbol <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

chrom_est <- get_firth_term(
  "GCAP_classChromosomal amplification"
)

ecdna_est <- get_firth_term(
  "GCAP_classecDNA"
)

tp53_est <- get_firth_term(
  "TP53_statusTP53 altered"
)

age_est <- get_firth_term(
  "Age_group>=10 years"
)

sex_est <- get_firth_term(
  "SexMale"
)

stage_est <- get_firth_term(
  "Stage_groupStage 3-4"
)

group_est <- get_firth_term(
  "Clinical_group_binaryGroup III-IV"
)

n_total_mv <- nrow(
  erms_multivariable_3group
)

n_no_focal <- sum(
  erms_multivariable_3group$GCAP_class ==
    "No focal amplification"
)

n_chrom <- sum(
  erms_multivariable_3group$GCAP_class ==
    "Chromosomal amplification"
)

n_ecdna <- sum(
  erms_multivariable_3group$GCAP_class ==
    "ecDNA"
)

n_tp53_wt <- sum(
  erms_multivariable_3group$TP53_status ==
    "TP53 WT"
)

n_tp53_alt <- sum(
  erms_multivariable_3group$TP53_status ==
    "TP53 altered"
)

n_female <- sum(
  erms_multivariable_3group$Sex ==
    "Female"
)

n_male <- sum(
  erms_multivariable_3group$Sex ==
    "Male"
)


n_stage12 <- sum(
  erms_multivariable_3group$Stage_group ==
    "Stage 1-2"
)

n_stage34 <- sum(
  erms_multivariable_3group$Stage_group ==
    "Stage 3-4"
)

n_group12 <- sum(
  erms_multivariable_3group$Clinical_group_binary ==
    "Group I-II"
)

n_group34 <- sum(
  erms_multivariable_3group$Clinical_group_binary ==
    "Group III-IV"
)

forest_header_row <- function(
    section,
    p_text = ""
) {
  
  tibble::tibble(
    Section = section,
    Level = "",
    Value_n = "",
    HR = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    P_value = NA_real_,
    HR_text = "",
    P_text = p_text,
    is_header = TRUE,
    is_reference = FALSE
  )
}

forest_model_row <- function(
    level,
    value_n,
    estimate,
    reference = FALSE
) {
  
  if (reference) {
    
    return(
      tibble::tibble(
        Section = "",
        Level = level,
        Value_n = value_n,
        HR = 1,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_,
        HR_text = "reference",
        P_text = "",
        is_header = FALSE,
        is_reference = TRUE
      )
    )
  }
  
  tibble::tibble(
    Section = "",
    Level = level,
    Value_n = value_n,
    HR = estimate$HR,
    CI_lower = estimate$CI_lower,
    CI_upper = estimate$CI_upper,
    P_value = estimate$P_value,
    HR_text = fmt_hr_ci(
      estimate$HR,
      estimate$CI_lower,
      estimate$CI_upper
    ),
    P_text = paste0(
      fmt_model_p(
        estimate$P_value
      ),
      " ",
      sig_symbol(
        estimate$P_value
      )
    ),
    is_header = FALSE,
    is_reference = FALSE
  )
}

forest_pub <- dplyr::bind_rows(
  
  forest_header_row(
    "Amplification class",
    paste0(
      "global P=",
      fmt_model_p(
        erms_3group_clinical_global_p
      )
    )
  ),
  
  forest_model_row(
    "No focal amplification",
    paste0("N=", n_no_focal),
    chrom_est,
    reference = TRUE
  ),
  
  forest_model_row(
    "Chromosomal amplification",
    paste0("N=", n_chrom),
    chrom_est
  ),
  
  forest_model_row(
    "ecDNA",
    paste0("N=", n_ecdna),
    ecdna_est
  ),
  
  forest_header_row(
    "TP53 status"
  ),
  
  forest_model_row(
    "WT",
    paste0("N=", n_tp53_wt),
    tp53_est,
    reference = TRUE
  ),
  
  forest_model_row(
    "Altered",
    paste0("N=", n_tp53_alt),
    tp53_est
  ),
  
  forest_header_row(
    "Sex"
  ),
  
  forest_model_row(
    "Female",
    paste0("N=", n_female),
    sex_est,
    reference = TRUE
  ),
  
  forest_model_row(
    "Male",
    paste0("N=", n_male),
    sex_est
  ),
  
  forest_header_row(
    "Age at diagnosis"
  ),
  
  forest_model_row(
    "<10 years",
    paste0(
      "N=",
      sum(
        erms_multivariable_3group$Age_group == "<10 years"
      )
    ),
    age_est,
    reference = TRUE
  ),
  
  forest_model_row(
    ">=10 years",
    paste0(
      "N=",
      sum(
        erms_multivariable_3group$Age_group == ">=10 years"
      )
    ),
    age_est
  ),
  
  forest_header_row(
    "Local stage"
  ),
  
  forest_model_row(
    "Stage 1-2",
    paste0("N=", n_stage12),
    stage_est,
    reference = TRUE
  ),
  
  forest_model_row(
    "Stage 3-4",
    paste0("N=", n_stage34),
    stage_est
  ),
  
  forest_header_row(
    "Clinical group"
  ),
  
  forest_model_row(
    "Group I-II",
    paste0("N=", n_group12),
    group_est,
    reference = TRUE
  ),
  
  forest_model_row(
    "Group III-IV",
    paste0("N=", n_group34),
    group_est
  )
)

forest_pub <- forest_pub |>
  dplyr::mutate(
    row_id = rev(
      seq_len(
        dplyr::n()
      )
    )
  )


header_rows <- forest_pub |>
  dplyr::filter(
    is_header
  )

# Wider table layout to avoid text overwriting.
table_xmin <- 0.50
table_xmax <- 5.35

x_parameter <- 0.62
x_level <- 0.86
x_value <- 2.70
x_hr <- 3.82
x_p <- 5.18

table_panel <- ggplot2::ggplot() +
  
  ggplot2::geom_rect(
    data = header_rows,
    ggplot2::aes(
      xmin = table_xmin,
      xmax = table_xmax,
      ymin = row_id - 0.43,
      ymax = row_id + 0.43
    ),
    fill = "grey90",
    color = NA
  ) +
  
  ggplot2::annotate(
    "text",
    x = c(
      x_parameter,
      x_value,
      x_hr,
      x_p
    ),
    y = max(
      forest_pub$row_id
    ) + 1.10,
    label = c(
      "parameter",
      "value (n)",
      "hazard ratio (95% CI)",
      "p"
    ),
    hjust = c(
      0,
      0,
      0.5,
      1
    ),
    fontface = "bold",
    size = 5.0
  ) +
  
  ggplot2::geom_text(
    data = forest_pub |>
      dplyr::filter(
        Section != ""
      ),
    ggplot2::aes(
      x = x_parameter,
      y = row_id,
      label = Section
    ),
    hjust = 0,
    fontface = "bold",
    size = 4.5
  ) +
  
  ggplot2::geom_text(
    data = forest_pub |>
      dplyr::filter(
        Level != ""
      ),
    ggplot2::aes(
      x = x_level,
      y = row_id,
      label = Level
    ),
    hjust = 0,
    size = 4.2
  ) +
  
  ggplot2::geom_text(
    data = forest_pub |>
      dplyr::filter(
        Value_n != ""
      ),
    ggplot2::aes(
      x = x_value,
      y = row_id,
      label = Value_n
    ),
    hjust = 0,
    size = 4.2
  ) +
  
  ggplot2::geom_text(
    data = forest_pub |>
      dplyr::filter(
        HR_text != ""
      ),
    ggplot2::aes(
      x = x_hr,
      y = row_id,
      label = HR_text
    ),
    hjust = 0.5,
    size = 4.2
  ) +
  
  ggplot2::geom_text(
    data = forest_pub |>
      dplyr::filter(
        P_text != ""
      ),
    ggplot2::aes(
      x = x_p,
      y = row_id,
      label = P_text
    ),
    hjust = 1,
    size = 4.2
  ) +
  
  ggplot2::coord_cartesian(
    xlim = c(
      table_xmin,
      table_xmax
    ),
    ylim = c(
      0.35,
      max(
        forest_pub$row_id
      ) + 1.50
    ),
    clip = "off"
  ) +
  
  ggplot2::theme_void(
    base_size = 13
  )

forest_estimates <- forest_pub |>
  dplyr::filter(
    !is.na(HR),
    HR > 0
  )

forest_ci <- forest_estimates |>
  dplyr::filter(
    !is_reference,
    !is.na(CI_lower),
    !is.na(CI_upper),
    CI_lower > 0,
    CI_upper > 0
  )

forest_xmin <- max(
  0.02,
  min(
    forest_ci$CI_lower,
    na.rm = TRUE
  ) / 1.5
)

forest_xmax <- max(
  10,
  max(
    forest_ci$CI_upper,
    na.rm = TRUE
  ) * 1.15
)

forest_panel <- ggplot2::ggplot() +
  
  ggplot2::geom_rect(
    data = header_rows,
    ggplot2::aes(
      xmin = forest_xmin,
      xmax = forest_xmax,
      ymin = row_id - 0.43,
      ymax = row_id + 0.43
    ),
    fill = "grey90",
    color = NA
  ) +
  
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.6
  ) +
  
  ggplot2::geom_segment(
    data = forest_ci,
    ggplot2::aes(
      x = CI_lower,
      xend = CI_upper,
      y = row_id,
      yend = row_id
    ),
    linewidth = 0.72
  ) +
  
  ggplot2::geom_point(
    data = forest_estimates |>
      dplyr::filter(
        !is_reference
      ),
    ggplot2::aes(
      x = HR,
      y = row_id
    ),
    shape = 15,
    size = 3.25
  ) +
  
  ggplot2::geom_point(
    data = forest_estimates |>
      dplyr::filter(
        is_reference
      ),
    ggplot2::aes(
      x = 1,
      y = row_id
    ),
    shape = 15,
    size = 2.9
  ) +
  
  ggplot2::scale_x_log10(
    breaks = c(
      0.05,
      0.1,
      0.25,
      0.5,
      1,
      2,
      5,
      10,
      20,
      50
    )
  ) +
  
  ggplot2::coord_cartesian(
    xlim = c(
      forest_xmin,
      forest_xmax
    ),
    ylim = c(
      0.35,
      max(
        forest_pub$row_id
      ) + 1.50
    ),
    clip = "off"
  ) +
  
  ggplot2::labs(
    x = "Adjusted hazard ratio",
    y = NULL
  ) +
  
  ggplot2::theme_classic(
    base_size = 13
  ) +
  
  ggplot2::theme(
    axis.text.y =
      ggplot2::element_blank(),
    axis.ticks.y =
      ggplot2::element_blank(),
    axis.line.y =
      ggplot2::element_blank(),
    panel.grid.major.x =
      ggplot2::element_line(
        linewidth = 0.35,
        color = "grey88"
      ),
    panel.grid.minor =
      ggplot2::element_blank()
  )

forest_grob <- gridExtra::arrangeGrob(
  table_panel,
  forest_panel,
  ncol = 2,
  widths = c(
    2.1,
    1
  ),
  top = grid::textGrob(
    paste0(
      "ERMS multivariable survival analysis   ",
      "Firth Cox regression; n=",
      nrow(
        erms_multivariable_3group
      ),
      "; deaths=",
      sum(
        erms_multivariable_3group$event,
        na.rm = TRUE
      )
    ),
    x = 0.01,
    hjust = 0,
    gp = grid::gpar(
      fontsize = 18,
      fontface = "bold"
    )
  )
)



grid::grid.newpage()
grid::grid.draw(
  forest_grob
)



grid::grid.newpage()

grid::grid.draw(
  forest_grob
)

cat(
  "\nERMS three-group model WITHOUT germline complete.\n"
)

cat(
  "Age is modeled as <10 years vs >=10 years.\n"
)


# =============================================================================
# 16. COMBINED ERMS KM + MULTIVARIABLE FOREST FIGURE
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Check that the forest plot has already been created
# -----------------------------------------------------------------------------

if (!exists("forest_grob")) {
  stop(
    "Object `forest_grob` was not found. ",
    "Run the ERMS multivariable/forest-plot section first."
  )
}

# -----------------------------------------------------------------------------
# 2. ERMS dataset for Kaplan-Meier analysis
# -----------------------------------------------------------------------------

erms_km_df <- gcap_patient |>
  dplyr::filter(
    Tumor_type == "ERMS",
    !is.na(GCAP_class),
    !is.na(Time_of_FU),
    !is.na(event)
  ) |>
  dplyr::mutate(
    GCAP_class = factor(
      GCAP_class,
      levels = c(
        "No focal amplification",
        "Chromosomal amplification",
        "ecDNA"
      )
    )
  ) |>
  droplevels()

cat(
  "\nERMS Kaplan-Meier cohort: n=",
  nrow(erms_km_df),
  "; deaths=",
  sum(erms_km_df$event, na.rm = TRUE),
  "\n",
  sep = ""
)

print(
  table(
    erms_km_df$GCAP_class,
    erms_km_df$event,
    useNA = "ifany"
  )
)

# -----------------------------------------------------------------------------
# 3. Kaplan-Meier fit + global log-rank P value
# -----------------------------------------------------------------------------

fit_erms_km <- survival::survfit(
  survival::Surv(
    Time_of_FU,
    event
  ) ~ GCAP_class,
  data = erms_km_df
)

logrank_erms <- survival::survdiff(
  survival::Surv(
    Time_of_FU,
    event
  ) ~ GCAP_class,
  data = erms_km_df
)

global_logrank_p <- stats::pchisq(
  logrank_erms$chisq,
  df = length(logrank_erms$n) - 1,
  lower.tail = FALSE
)

format_p_combined <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  
  if (p < 0.001) {
    return("<0.001")
  }
  
  sprintf("%.3f", p)
}

cat(
  "Global ERMS log-rank P = ",
  format_p_combined(global_logrank_p),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# 4. Colors
# -----------------------------------------------------------------------------

km_palette <- c(
  "No focal amplification" = "#008B45FF",
  "Chromosomal amplification" = "#EE0000FF",
  "ecDNA" = "#3B4992FF"
)

# -----------------------------------------------------------------------------
# 5. Publication-style KM panel
# -----------------------------------------------------------------------------

km_plot_object <- survminer::ggsurvplot(
  fit_erms_km,
  data = erms_km_df,
  
  conf.int = FALSE,
  
  pval = paste0(
    "Global log-rank P=",
    format_p_combined(global_logrank_p)
  ),
  
  pval.coord = c(3, 0.20),
  pval.size = 5.5,
  
  risk.table = TRUE,
  risk.table.height = 0.25,
  
  censor = TRUE,
  censor.shape = 124,
  censor.size = 2.5,
  
  xlab = "Follow-up time (months)",
  ylab = "Overall survival",
  
  xlim = c(
    0,
    max(erms_km_df$Time_of_FU, na.rm = TRUE)
  ),
  
  break.time.by = 40,
  ylim = c(0, 1.05),
  
  palette = unname(
    km_palette[
      levels(erms_km_df$GCAP_class)
    ]
  ),
  
  legend.title = NULL,
  legend.labs = levels(erms_km_df$GCAP_class),
  legend = "bottom",
  
  ggtheme = ggplot2::theme_classic(
    base_size = 14
  )
)

# Main KM plot
km_main <- km_plot_object$plot +
  ggplot2::labs(tag = "a") +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 16
    ),
    plot.tag.position = c(0.01, 0.99),
    axis.title = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 11.5),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 8, 2, 8)
  )

# Risk table
km_risk_table <- km_plot_object$table +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 10),
    legend.position = "none",
    plot.margin = ggplot2::margin(0, 8, 8, 8)
  )

# Combine KM curve + risk table into one panel
km_panel <- gridExtra::arrangeGrob(
  ggplot2::ggplotGrob(km_main),
  ggplot2::ggplotGrob(km_risk_table),
  ncol = 1,
  heights = c(0.76, 0.24)
)

# -----------------------------------------------------------------------------
# 6. Add panel label "b" to the forest plot
#
# Because forest_grob is already a grob, put it inside another grob and add
# the panel label separately.
# -----------------------------------------------------------------------------

forest_panel_with_tag <- gridExtra::arrangeGrob(
  grid::textGrob(
    "b",
    x = 0,
    y = 1,
    just = c(
      "left",
      "top"
    ),
    gp = grid::gpar(
      fontsize = 18,
      fontface = "bold"
    )
  ),
  forest_grob,
  
  ncol = 1,
  
  heights = c(
    0.045,
    0.955
  )
)

# -----------------------------------------------------------------------------
# 7. Combine KM + forest
# -----------------------------------------------------------------------------

combined_erms_figure <- gridExtra::arrangeGrob(
  km_panel,
  forest_panel_with_tag,
  
  ncol = 2,
  
  # Give the forest plot more horizontal space.
  widths = c(
    0.95,
    2.15
  ),
  
  top = grid::textGrob(
    "Overall survival of ERMS patients by amplification architecture",
    
    x = 0.01,
    hjust = 0,
    
    gp = grid::gpar(
      fontsize = 17,
      fontface = "bold"
    )
  )
)

# -----------------------------------------------------------------------------
# 8. Display combined figure
# -----------------------------------------------------------------------------



grid::grid.newpage()

grid::grid.draw(
  combined_erms_figure
)



# Display in RStudio.
grid::grid.newpage()

grid::grid.draw(
  combined_erms_figure
)


# =============================================================================
# 17. RECURRENT CYTOBANDS AND IMPACT GENES BY AMPLICON CLASS
#
# Uses the GCAP-filtered sample set already created above (`gcap`).
#
# Analyses:
#   1. Recurrent circular/ecDNA cytobands — sample level
#   2. Recurrent circular/ecDNA cytobands — patient level
#   3. Recurrent IMPACT genes on circular/ecDNA — sample level
#   4. Recurrent IMPACT genes on circular/ecDNA — patient level
#   5. Circular vs noncircular amplified IMPACT genes — sample level
#   6. Circular vs noncircular amplified IMPACT genes — patient level
#
# Recurrence threshold:
#   >= 2 independent samples or patients.
#
# NOTE:
# =============================================================================



min_recurrent_samples <- 2L
min_recurrent_patients <- 2L

# Match the amplification-architecture colors used elsewhere in this script.
amplicon_class_colors <- c(
  "noncircular" = "#EE0000FF",
  "circular" = "#3B4992FF"
)

amplicon_class_labels <- c(
  "noncircular" = "Chromosomal amplification",
  "circular" = "ecDNA"
)

# -----------------------------------------------------------------------------
# 17A. Helpers
# -----------------------------------------------------------------------------

require_columns <- function(
    data,
    required,
    data_name
) {
  
  missing_cols <- setdiff(
    required,
    names(data)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      data_name,
      " is missing required column(s): ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }
}

collapse_nonmissing <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x <- sort(
    unique(
      x[
        !is.na(x) &
          x != ""
      ]
    )
  )
  
  if (length(x) == 0) {
    return("")
  }
  
  paste(
    x,
    collapse = "; "
  )
}

first_nonmissing_recurrent <- function(x) {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      trimws(x) != ""
  ]
  
  if (length(x) == 0) {
    return(
      NA_character_
    )
  }
  
  x[1]
}

integer_breaks <- function(x) {
  
  xmax <- suppressWarnings(
    max(
      x,
      na.rm = TRUE
    )
  )
  
  if (
    !is.finite(xmax) ||
    xmax < 1
  ) {
    xmax <- 1
  }
  
  seq(
    0,
    ceiling(xmax),
    by = 1
  )
}


# -----------------------------------------------------------------------------
# 17B. Read and prepare recurrent-gene data
# -----------------------------------------------------------------------------

genes_gcap <- read.csv(
  recurrent_genes_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

impact_genes <- read.csv(
  impact_gene_list_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# dplyr cannot operate on data frames containing blank or NA column names.
# Repair names immediately after import.
names(genes_gcap) <- make.names(
  names(genes_gcap),
  unique = TRUE
)

names(impact_genes) <- make.names(
  names(impact_genes),
  unique = TRUE
)

require_columns(
  genes_gcap,
  c(
    "sample",
    "gene_class",
    "band",
    "hgnc_symbol"
  ),
  "GCAP annotated gene file"
)

require_columns(
  impact_genes,
  "Gene",
  "IMPACT gene list"
)

require_columns(
  gcap,
  c(
    "sample",
    "Patient.ID",
    "Tumor_type"
  ),
  "GCAP sample table"
)

target_genes <- impact_genes |>
  dplyr::transmute(
    Gene = trimws(
      as.character(Gene)
    )
  ) |>
  dplyr::filter(
    !is.na(Gene),
    Gene != ""
  ) |>
  dplyr::distinct(
    Gene
  ) |>
  dplyr::pull(
    Gene
  )

# One metadata row per GCAP-filtered sample.
gcap_sample_metadata <- gcap |>
  dplyr::transmute(
    sample = as.character(sample),
    Patient.ID = as.character(Patient.ID),
    Tumor_type = as.character(Tumor_type)
  ) |>
  dplyr::filter(
    !is.na(sample),
    sample != ""
  ) |>
  dplyr::group_by(
    sample
  ) |>
  dplyr::summarise(
    Patient.ID = first_nonmissing_recurrent(
      Patient.ID
    ),
    Tumor_type = first_nonmissing_recurrent(
      Tumor_type
    ),
    .groups = "drop"
  )

genes_gcap_sub <- genes_gcap |>
  dplyr::mutate(
    sample = as.character(sample),
    gene_class = tolower(
      trimws(
        as.character(gene_class)
      )
    ),
    band = trimws(
      as.character(band)
    ),
    hgnc_symbol = trimws(
      as.character(hgnc_symbol)
    )
  ) |>
  dplyr::inner_join(
    gcap_sample_metadata,
    by = "sample"
  )

genes_circular <- genes_gcap_sub |>
  dplyr::filter(
    gene_class == "circular"
  )

genes_amplified <- genes_gcap_sub |>
  dplyr::filter(
    gene_class %in% c(
      "circular",
      "noncircular"
    )
  )

cat(
  "\nRecurrent-gene analysis:\n",
  "  GCAP-filtered samples represented = ",
  dplyr::n_distinct(
    genes_gcap_sub$sample
  ),
  "\n",
  "  Patients represented = ",
  dplyr::n_distinct(
    genes_gcap_sub$Patient.ID,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

# =============================================================================
# 17C. RECURRENT CIRCULAR/ecDNA CYTOBANDS — SAMPLE LEVEL
# =============================================================================

circular_band_sample <- genes_circular |>
  dplyr::filter(
    !is.na(band),
    band != ""
  ) |>
  dplyr::distinct(
    sample,
    band
  ) |>
  dplyr::count(
    band,
    name = "n_samples"
  ) |>
  dplyr::filter(
    n_samples >=
      min_recurrent_samples
  ) |>
  dplyr::arrange(
    dplyr::desc(
      n_samples
    ),
    band
  ) |>
  dplyr::mutate(
    band = stats::reorder(
      band,
      n_samples
    )
  )

if (nrow(circular_band_sample) > 0) {
  
  p_circular_band_sample <- ggplot2::ggplot(
    circular_band_sample,
    ggplot2::aes(
      x = band,
      y = n_samples
    )
  ) +
    ggplot2::geom_col(
      fill = "#6C77AD"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        circular_band_sample$n_samples
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::labs(
      title =
        "Recurrent cytobands on ecDNA — sample level",
      x = "Cytoband",
      y = "Number of samples"
    )
  
  print(
    p_circular_band_sample
  )
  
}


# =============================================================================
# 17D. RECURRENT CIRCULAR/ecDNA CYTOBANDS — PATIENT LEVEL
# =============================================================================

circular_band_patient <- genes_circular |>
  dplyr::filter(
    !is.na(Patient.ID),
    Patient.ID != "",
    !is.na(band),
    band != ""
  ) |>
  dplyr::distinct(
    Patient.ID,
    band
  ) |>
  dplyr::count(
    band,
    name = "n_patients"
  ) |>
  dplyr::filter(
    n_patients >=
      min_recurrent_patients
  ) |>
  dplyr::arrange(
    dplyr::desc(
      n_patients
    ),
    band
  ) |>
  dplyr::mutate(
    band = stats::reorder(
      band,
      n_patients
    )
  )

if (nrow(circular_band_patient) > 0) {
  
  p_circular_band_patient <- ggplot2::ggplot(
    circular_band_patient,
    ggplot2::aes(
      x = band,
      y = n_patients
    )
  ) +
    ggplot2::geom_col(
      fill = "#6C77AD"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        circular_band_patient$n_patients
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::labs(
      title =
        "Recurrent cytobands on ecDNA — patient level",
      x = "Cytoband",
      y = "Number of patients"
    )
  
  print(
    p_circular_band_patient
  )
  
}


# =============================================================================
# 17E. RECURRENT IMPACT GENES ON CIRCULAR/ecDNA — SAMPLE LEVEL
# =============================================================================

circular_gene_sample <- genes_circular |>
  dplyr::filter(
    hgnc_symbol %in%
      target_genes,
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) |>
  dplyr::group_by(
    hgnc_symbol
  ) |>
  dplyr::summarise(
    band = collapse_nonmissing(
      band
    ),
    n_samples = dplyr::n_distinct(
      sample
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_samples >=
      min_recurrent_samples
  ) |>
  dplyr::arrange(
    dplyr::desc(
      n_samples
    ),
    hgnc_symbol
  ) |>
  dplyr::mutate(
    gene_band = dplyr::if_else(
      band == "",
      hgnc_symbol,
      paste0(
        hgnc_symbol,
        " (",
        band,
        ")"
      )
    ),
    gene_band = stats::reorder(
      gene_band,
      n_samples
    )
  )

if (nrow(circular_gene_sample) > 0) {
  
  p_circular_gene_sample <- ggplot2::ggplot(
    circular_gene_sample,
    ggplot2::aes(
      x = gene_band,
      y = n_samples
    )
  ) +
    ggplot2::geom_col(
      fill = "#6C77AD"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        circular_gene_sample$n_samples
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::labs(
      title =
        "Recurrent IMPACT genes on ecDNA — sample level",
      x = "Gene (cytoband)",
      y = "Number of samples"
    )
  
  print(
    p_circular_gene_sample
  )
  
}


# =============================================================================
# 17F. RECURRENT IMPACT GENES ON CIRCULAR/ecDNA — PATIENT LEVEL
# =============================================================================
amplicon_class_colors <- c(
  "noncircular" = "#EE0000",
  "circular" = "#6C77AD"
)

amplicon_class_labels <- c(
  "noncircular" = "Chromosomal",
  "circular" = "ecDNA"
)
circular_gene_patient <- genes_circular |>
  dplyr::filter(
    hgnc_symbol %in%
      target_genes,
    !is.na(hgnc_symbol),
    hgnc_symbol != "",
    !is.na(Patient.ID),
    Patient.ID != ""
  ) |>
  dplyr::group_by(
    hgnc_symbol
  ) |>
  dplyr::summarise(
    band = collapse_nonmissing(
      band
    ),
    n_patients = dplyr::n_distinct(
      Patient.ID
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_patients >=
      min_recurrent_patients
  ) |>
  dplyr::arrange(
    dplyr::desc(
      n_patients
    ),
    hgnc_symbol
  ) |>
  dplyr::mutate(
    gene_band = dplyr::if_else(
      band == "",
      hgnc_symbol,
      paste0(
        hgnc_symbol,
        " (",
        band,
        ")"
      )
    ),
    gene_band = stats::reorder(
      gene_band,
      n_patients
    )
  )

if (nrow(circular_gene_patient) > 0) {
  
  p_circular_gene_patient <- ggplot2::ggplot(
    circular_gene_patient,
    ggplot2::aes(
      x = gene_band,
      y = n_patients
    )
  ) +
    ggplot2::geom_col(
      fill = amplicon_class_colors[["circular"]]
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        circular_gene_patient$n_patients
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::labs(
      title =
        "Recurrent IMPACT genes on ecDNA — patient level",
      x = "Gene (cytoband)",
      y = "Number of patients"
    )
  
  print(
    p_circular_gene_patient
  )
  
}


# =============================================================================
# 17G. CIRCULAR VS NONCIRCULAR IMPACT GENES — SAMPLE LEVEL
# =============================================================================

amplified_impact_genes <- genes_amplified |>
  dplyr::filter(
    hgnc_symbol %in%
      target_genes,
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  )

# Select recurrent genes based on UNIQUE samples carrying the amplified gene,
# independent of architecture. This prevents double counting when choosing
# which genes meet the recurrence threshold.
top_genes_sample <- amplified_impact_genes |>
  dplyr::distinct(
    sample,
    hgnc_symbol
  ) |>
  dplyr::count(
    hgnc_symbol,
    name = "total_samples"
  ) |>
  dplyr::filter(
    total_samples >=
      min_recurrent_samples
  ) |>
  dplyr::arrange(
    dplyr::desc(
      total_samples
    ),
    hgnc_symbol
  )

gene_band_lookup <- amplified_impact_genes |>
  dplyr::group_by(
    hgnc_symbol
  ) |>
  dplyr::summarise(
    band = collapse_nonmissing(
      band
    ),
    .groups = "drop"
  )

gene_class_sample_freq <- amplified_impact_genes |>
  dplyr::filter(
    hgnc_symbol %in%
      top_genes_sample$hgnc_symbol
  ) |>
  dplyr::distinct(
    sample,
    hgnc_symbol,
    gene_class
  ) |>
  dplyr::count(
    hgnc_symbol,
    gene_class,
    name = "n_samples"
  ) |>
  dplyr::left_join(
    top_genes_sample,
    by = "hgnc_symbol"
  ) |>
  dplyr::left_join(
    gene_band_lookup,
    by = "hgnc_symbol"
  ) |>
  dplyr::mutate(
    gene_band = dplyr::if_else(
      band == "",
      hgnc_symbol,
      paste0(
        hgnc_symbol,
        " (",
        band,
        ")"
      )
    ),
    gene_band = stats::reorder(
      gene_band,
      total_samples
    ),
    gene_class = factor(
      gene_class,
      levels = c(
        "noncircular",
        "circular"
      )
    )
  )

if (nrow(gene_class_sample_freq) > 0) {
  
  p_gene_class_sample <- ggplot2::ggplot(
    gene_class_sample_freq,
    ggplot2::aes(
      x = gene_band,
      y = n_samples,
      fill = gene_class
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values =
        amplicon_class_colors,
      labels =
        amplicon_class_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        gene_class_sample_freq$n_samples
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::labs(
      title =
        "Recurrent amplified IMPACT genes by architecture — sample level",
      x = "Gene (cytoband)",
      y = "Number of samples",
      fill = "Amplicon class"
    )
  
  print(
    p_gene_class_sample
  )
  
}


# =============================================================================
# 17H. CIRCULAR VS NONCIRCULAR IMPACT GENES — PATIENT LEVEL
# =============================================================================

top_genes_patient <- amplified_impact_genes |>
  dplyr::filter(
    !is.na(Patient.ID),
    Patient.ID != ""
  ) |>
  dplyr::distinct(
    Patient.ID,
    hgnc_symbol
  ) |>
  dplyr::count(
    hgnc_symbol,
    name = "total_patients"
  ) |>
  dplyr::filter(
    total_patients >=
      min_recurrent_patients
  ) |>
  dplyr::arrange(
    dplyr::desc(
      total_patients
    ),
    hgnc_symbol
  )

gene_class_patient_freq <- amplified_impact_genes |>
  dplyr::filter(
    !is.na(Patient.ID),
    Patient.ID != "",
    hgnc_symbol %in%
      top_genes_patient$hgnc_symbol
  ) |>
  dplyr::distinct(
    Patient.ID,
    hgnc_symbol,
    gene_class
  ) |>
  dplyr::count(
    hgnc_symbol,
    gene_class,
    name = "n_patients"
  ) |>
  dplyr::left_join(
    top_genes_patient,
    by = "hgnc_symbol"
  ) |>
  dplyr::left_join(
    gene_band_lookup,
    by = "hgnc_symbol"
  ) |>
  dplyr::mutate(
    gene_band = dplyr::if_else(
      band == "",
      hgnc_symbol,
      paste0(
        hgnc_symbol,
        " (",
        band,
        ")"
      )
    ),
    gene_band = stats::reorder(
      gene_band,
      total_patients
    ),
    gene_class = factor(
      gene_class,
      levels = c(
        "noncircular",
        "circular"
      )
    )
  )

if (nrow(gene_class_patient_freq) > 0) {
  
  p_gene_class_patient <- ggplot2::ggplot(
    gene_class_patient_freq,
    ggplot2::aes(
      x = gene_band,
      y = n_patients,
      fill = gene_class
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values =
        amplicon_class_colors,
      labels =
        amplicon_class_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        gene_class_patient_freq$n_patients
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::labs(
      title =
        "Recurrent amplified IMPACT genes by architecture — patient level",
      x = "Gene (cytoband)",
      y = "Number of patients",
      fill = "Amplicon class"
    )
  
  print(
    p_gene_class_patient
  )
  
}



# =============================================================================
# 18. RECURRENT ecDNA IMPACT GENES BY RMS SUBTYPE
#
# Genes are selected if they occur on circular/ecDNA amplification in at least:
#   - 2 independent patients (patient-level plot), or
#   - 2 independent samples  (sample-level plot).
#
# For those selected genes, BOTH chromosomal/noncircular and circular/ecDNA
# amplification counts are displayed within each RMS subtype.
#
# Uses objects created in Section 17:
#   genes_amplified
#   target_genes
#   amplicon_class_colors
#   amplicon_class_labels
#   collapse_nonmissing()
# =============================================================================

# -----------------------------------------------------------------------------
# 18A. Common gene/cytoband lookup
# -----------------------------------------------------------------------------

impact_gene_band_lookup <- genes_amplified |>
  dplyr::filter(
    hgnc_symbol %in% target_genes,
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) |>
  dplyr::group_by(
    hgnc_symbol
  ) |>
  dplyr::summarise(
    band = collapse_nonmissing(
      band
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    gene_band = dplyr::if_else(
      band == "",
      hgnc_symbol,
      paste0(
        hgnc_symbol,
        " (",
        band,
        ")"
      )
    )
  )

# =============================================================================
# 18B. PATIENT LEVEL
# =============================================================================

# Genes present on ecDNA in >=2 independent patients across the cohort.
circular_genes_patient <- genes_amplified |>
  dplyr::filter(
    gene_class == "circular",
    hgnc_symbol %in% target_genes,
    !is.na(Patient.ID),
    Patient.ID != "",
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) |>
  dplyr::distinct(
    Patient.ID,
    hgnc_symbol
  ) |>
  dplyr::count(
    hgnc_symbol,
    name = "n_circular_patients"
  ) |>
  dplyr::filter(
    n_circular_patients >=
      min_recurrent_patients
  )

# Patient counts by subtype, gene, and architecture.
plot_patient_type <- genes_amplified |>
  dplyr::filter(
    hgnc_symbol %in%
      circular_genes_patient$hgnc_symbol,
    gene_class %in%
      c(
        "noncircular",
        "circular"
      ),
    !is.na(Patient.ID),
    Patient.ID != "",
    !is.na(Tumor_type),
    Tumor_type != ""
  ) |>
  dplyr::distinct(
    Patient.ID,
    Tumor_type,
    hgnc_symbol,
    gene_class
  ) |>
  dplyr::count(
    Tumor_type,
    hgnc_symbol,
    gene_class,
    name = "n_patients"
  ) |>
  tidyr::complete(
    Tumor_type,
    hgnc_symbol =
      circular_genes_patient$hgnc_symbol,
    gene_class =
      c(
        "noncircular",
        "circular"
      ),
    fill = list(
      n_patients = 0
    )
  ) |>
  dplyr::left_join(
    circular_genes_patient,
    by = "hgnc_symbol"
  ) |>
  dplyr::left_join(
    impact_gene_band_lookup,
    by = "hgnc_symbol"
  ) |>
  dplyr::mutate(
    gene_class = factor(
      gene_class,
      levels = c(
        "noncircular",
        "circular"
      )
    ),
    gene_band = factor(
      gene_band,
      levels =
        impact_gene_band_lookup |>
        dplyr::filter(
          hgnc_symbol %in%
            circular_genes_patient$hgnc_symbol
        ) |>
        dplyr::left_join(
          circular_genes_patient,
          by = "hgnc_symbol"
        ) |>
        dplyr::arrange(
          n_circular_patients,
          hgnc_symbol
        ) |>
        dplyr::pull(
          gene_band
        ) |>
        unique()
    )
  )

if (nrow(plot_patient_type) > 0) {
  
  p_patient_type <- ggplot2::ggplot(
    plot_patient_type,
    ggplot2::aes(
      x = gene_band,
      y = n_patients,
      fill = gene_class
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(
      ~ Tumor_type,
      scales = "fixed"
    ) +
    ggplot2::scale_fill_manual(
      values =
        amplicon_class_colors,
      labels =
        amplicon_class_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        plot_patient_type$n_patients
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        face = "bold"
      )
    ) +
    ggplot2::labs(
      title =
        "Recurrent amplified IMPACT genes by RMS subtype — patient level",
      subtitle =
        "Genes shown have ecDNA amplification in at least two patients",
      x =
        "Gene (cytoband)",
      y =
        "Number of patients",
      fill =
        "Amplicon class"
    )
  
  print(
    p_patient_type
  )
  
}



# =============================================================================
# 18C. SAMPLE LEVEL
# =============================================================================

# Genes present on ecDNA in >=2 independent samples across the cohort.
circular_genes_sample <- genes_amplified |>
  dplyr::filter(
    gene_class == "circular",
    hgnc_symbol %in% target_genes,
    !is.na(sample),
    sample != "",
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) |>
  dplyr::distinct(
    sample,
    hgnc_symbol
  ) |>
  dplyr::count(
    hgnc_symbol,
    name = "n_circular_samples"
  ) |>
  dplyr::filter(
    n_circular_samples >=
      min_recurrent_samples
  )

# Sample counts by subtype, gene, and architecture.
plot_sample_type <- genes_amplified |>
  dplyr::filter(
    hgnc_symbol %in%
      circular_genes_sample$hgnc_symbol,
    gene_class %in%
      c(
        "noncircular",
        "circular"
      ),
    !is.na(sample),
    sample != "",
    !is.na(Tumor_type),
    Tumor_type != ""
  ) |>
  dplyr::distinct(
    sample,
    Tumor_type,
    hgnc_symbol,
    gene_class
  ) |>
  dplyr::count(
    Tumor_type,
    hgnc_symbol,
    gene_class,
    name = "n_samples"
  ) |>
  tidyr::complete(
    Tumor_type,
    hgnc_symbol =
      circular_genes_sample$hgnc_symbol,
    gene_class =
      c(
        "noncircular",
        "circular"
      ),
    fill = list(
      n_samples = 0
    )
  ) |>
  dplyr::left_join(
    circular_genes_sample,
    by = "hgnc_symbol"
  ) |>
  dplyr::left_join(
    impact_gene_band_lookup,
    by = "hgnc_symbol"
  ) |>
  dplyr::mutate(
    gene_class = factor(
      gene_class,
      levels = c(
        "noncircular",
        "circular"
      )
    ),
    gene_band = factor(
      gene_band,
      levels =
        impact_gene_band_lookup |>
        dplyr::filter(
          hgnc_symbol %in%
            circular_genes_sample$hgnc_symbol
        ) |>
        dplyr::left_join(
          circular_genes_sample,
          by = "hgnc_symbol"
        ) |>
        dplyr::arrange(
          n_circular_samples,
          hgnc_symbol
        ) |>
        dplyr::pull(
          gene_band
        ) |>
        unique()
    )
  )

if (nrow(plot_sample_type) > 0) {
  
  p_sample_type <- ggplot2::ggplot(
    plot_sample_type,
    ggplot2::aes(
      x = gene_band,
      y = n_samples,
      fill = gene_class
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(
      ~ Tumor_type,
      scales = "fixed"
    ) +
    ggplot2::scale_fill_manual(
      values =
        amplicon_class_colors,
      labels =
        amplicon_class_labels,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = integer_breaks(
        plot_sample_type$n_samples
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
    ) +
    ggplot2::theme_classic(
      base_size = 14
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        face = "bold"
      )
    ) +
    ggplot2::labs(
      title =
        "Recurrent amplified IMPACT genes by RMS subtype — sample level",
      subtitle =
        "Genes shown have ecDNA amplification in at least two samples",
      x =
        "Gene (cytoband)",
      y =
        "Number of samples",
      fill =
        "Amplicon class"
    )
  
  print(
    p_sample_type
  )
  
}



cat(
  "\nSubtype-stratified recurrent IMPACT-gene plots complete.\n",
  sep = ""
)
# =============================================================================
# MSK-IMPACT COPY-NUMBER FOLD CHANGE
# ALL AMPLIFIED GENES: ecDNA/circular vs chromosomal/noncircular
# =============================================================================

library(dplyr)
library(ggplot2)

# Load MSK-IMPACT fold-change table
fc <- read.csv(
  fold_change_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Consistent manuscript colors
fc_class_colors <- c(
  "noncircular" = "#EE0000FF",
  "circular" = "#6C77AD"
)

fc_class_labels <- c(
  "noncircular" = "Chromosomal",
  "circular" = "ecDNA"
)


# -----------------------------------------------------------------------------
# 1. Collapse GCAP architecture to one class per sample-gene
#    Uses ALL amplified genes
#    Priority: circular/ecDNA > noncircular
# -----------------------------------------------------------------------------
genes_amp <- genes_gcap_sub %>%
  filter(gene_class %in% c("circular", "noncircular"))

gene_class_df <- genes_amp |>
  dplyr::filter(
    !is.na(sample),
    sample != "",
    !is.na(hgnc_symbol),
    hgnc_symbol != "",
    gene_class %in% c(
      "circular",
      "noncircular"
    )
  ) |>
  dplyr::group_by(
    sample,
    hgnc_symbol
  ) |>
  dplyr::summarise(
    Patient.ID = dplyr::first(
      stats::na.omit(Patient.ID),
      default = NA_character_
    ),
    
    Tumor_type = dplyr::first(
      stats::na.omit(Tumor_type),
      default = NA_character_
    ),
    
    amp_class = dplyr::case_when(
      any(
        gene_class == "circular",
        na.rm = TRUE
      ) ~ "circular",
      
      any(
        gene_class == "noncircular",
        na.rm = TRUE
      ) ~ "noncircular",
      
      TRUE ~ NA_character_
    ),
    
    .groups = "drop"
  ) |>
  dplyr::filter(
    !is.na(amp_class)
  )


# -----------------------------------------------------------------------------
# 2. Clean MSK-IMPACT FoldChange table
#    Uses ALL genes called as amplification
# -----------------------------------------------------------------------------

impact_fc_clean <- fc |>
  dplyr::filter(
    Type == "Amplification"
  ) |>
  dplyr::transmute(
    sample = DMP_SAMPLE_ID,
    hgnc_symbol = Gene,
    FoldChange = as.numeric(
      FoldChange
    )
  ) |>
  dplyr::filter(
    !is.na(sample),
    sample != "",
    !is.na(hgnc_symbol),
    hgnc_symbol != "",
    !is.na(FoldChange),
    is.finite(FoldChange),
    FoldChange > 0
  ) |>
  dplyr::distinct()


# -----------------------------------------------------------------------------
# 3. Match FoldChange to amplification architecture
# -----------------------------------------------------------------------------

fc_gene_class <- impact_fc_clean |>
  dplyr::inner_join(
    gene_class_df,
    by = c(
      "sample",
      "hgnc_symbol"
    )
  ) |>
  dplyr::mutate(
    amp_class = factor(
      amp_class,
      levels = c(
        "noncircular",
        "circular"
      )
    ),
    log2_FC = log2(
      FoldChange
    )
  )


# -----------------------------------------------------------------------------
# 4. Summary counts
# -----------------------------------------------------------------------------

fc_counts <- fc_gene_class |>
  dplyr::summarise(
    n_gene_sample_observations = dplyr::n(),
    n_samples = dplyr::n_distinct(
      sample
    ),
    n_patients = dplyr::n_distinct(
      Patient.ID
    ),
    n_genes = dplyr::n_distinct(
      hgnc_symbol
    )
  )

fc_counts_by_class <- fc_gene_class |>
  dplyr::group_by(
    amp_class
  ) |>
  dplyr::summarise(
    n_gene_sample_observations = dplyr::n(),
    n_samples = dplyr::n_distinct(
      sample
    ),
    n_patients = dplyr::n_distinct(
      Patient.ID
    ),
    n_genes = dplyr::n_distinct(
      hgnc_symbol
    ),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 5. Descriptive statistics
# -----------------------------------------------------------------------------

fc_summary <- fc_gene_class |>
  dplyr::group_by(
    amp_class
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    
    median_FC = median(
      FoldChange,
      na.rm = TRUE
    ),
    
    Q1_FC = stats::quantile(
      FoldChange,
      0.25,
      na.rm = TRUE
    ),
    
    Q3_FC = stats::quantile(
      FoldChange,
      0.75,
      na.rm = TRUE
    ),
    
    mean_FC = mean(
      FoldChange,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 6. Wilcoxon rank-sum test
# -----------------------------------------------------------------------------

fc_wilcox <- stats::wilcox.test(
  FoldChange ~ amp_class,
  data = fc_gene_class,
  exact = FALSE
)

fc_p <- fc_wilcox$p.value

fc_p_label <- ifelse(
  fc_p < 0.001,
  "P < 0.001",
  paste0(
    "P = ",
    sprintf(
      "%.3f",
      fc_p
    )
  )
)


# -----------------------------------------------------------------------------
# 7. Box plot
# -----------------------------------------------------------------------------

fc_ymax <- max(
  fc_gene_class$FoldChange,
  na.rm = TRUE
)

p_fc_overall <- ggplot2::ggplot(
  fc_gene_class,
  ggplot2::aes(
    x = amp_class,
    y = FoldChange,
    fill = amp_class
  )
) +
  ggplot2::geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(
      color = amp_class
    ),
    width = 0.12,
    size = 1.8,
    alpha = 0.7
  ) +
  ggplot2::annotate(
    "segment",
    x = 1,
    xend = 2,
    y = fc_ymax * 1.08,
    yend = fc_ymax * 1.08,
    linewidth = 0.4
  ) +
  ggplot2::annotate(
    "segment",
    x = 1,
    xend = 1,
    y = fc_ymax * 1.05,
    yend = fc_ymax * 1.08,
    linewidth = 0.4
  ) +
  ggplot2::annotate(
    "segment",
    x = 2,
    xend = 2,
    y = fc_ymax * 1.05,
    yend = fc_ymax * 1.08,
    linewidth = 0.4
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = fc_ymax * 1.13,
    label = fc_p_label,
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = fc_class_colors,
    labels = fc_class_labels
  ) +
  ggplot2::scale_color_manual(
    values = fc_class_colors,
    labels = fc_class_labels
  ) +
  ggplot2::scale_x_discrete(
    labels = fc_class_labels
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.05,
        0.20
      )
    )
  ) +
  ggplot2::theme_classic(
    base_size = 14
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 12
    )
  ) +
  ggplot2::labs(
    y = "MSK-IMPACT copy-number fold change"
  )

fc_counts
fc_counts_by_class
fc_summary
fc_wilcox
p_fc_overall





# =============================================================================
# 19. 12q13-15 ARCHITECTURE SURVIVAL
#
# Driver-locus definition:
#   MDM2, CDK4, GLI1
#
# Patient-level hierarchy:
#   12q ecDNA > Chromosomal 12q amplification > No focal 12q amplification
#
# Univariable analyses:
#   - Overall Firth Cox adjusted for RMS subtype
#   - ARMS / ERMS / SSRMS subtype-specific Firth Cox
#
# ERMS multivariable analysis:
#   12q architecture + TP53 status + age + sex + local stage + clinical group
#
# This section REUSES erms_multivariable_3group_all from Section 15.
# It does not rebuild the clinical covariates and does not run TP53-only models.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(coxphf)
  library(ggplot2)
  library(cowplot)
})

# =============================================================================
# 1. DEFINE 12q13-15 DRIVER LOCUS
# =============================================================================

q12_driver_genes <- c(
  "MDM2",
  "CDK4",
  "GLI1"
)

q12_architecture_levels <- c(
  "No focal 12q amplification",
  "Chromosomal 12q amplification",
  "12q ecDNA"
)

q12_architecture_colors <- c(
  "No focal 12q amplification" = "#008B45FF",
  "Chromosomal 12q amplification" = "#EE0000FF",
  "12q ecDNA" = "#3B4992FF"
)

# =============================================================================
# 2. BUILD PATIENT-LEVEL 12q13-15 ARCHITECTURE
# =============================================================================

q12_driver_records <- genes_amplified |>
  dplyr::filter(
    hgnc_symbol %in% q12_driver_genes,
    !is.na(Patient.ID),
    Patient.ID != "",
    gene_class %in% c(
      "circular",
      "noncircular"
    )
  ) |>
  dplyr::distinct(
    Patient.ID,
    sample,
    Tumor_type,
    hgnc_symbol,
    band,
    gene_class
  )

q12_patient_observed <- q12_driver_records |>
  dplyr::group_by(
    Patient.ID
  ) |>
  dplyr::summarise(
    q12_any_amplification = 1L,
    
    q12_any_ecDNA = as.integer(
      any(
        gene_class == "circular",
        na.rm = TRUE
      )
    ),
    
    q12_any_chromosomal = as.integer(
      any(
        gene_class == "noncircular",
        na.rm = TRUE
      )
    ),
    
    q12_genes = paste(
      sort(
        unique(
          hgnc_symbol[
            !is.na(hgnc_symbol) &
              hgnc_symbol != ""
          ]
        )
      ),
      collapse = "; "
    ),
    
    q12_bands = paste(
      sort(
        unique(
          band[
            !is.na(band) &
              band != ""
          ]
        )
      ),
      collapse = "; "
    ),
    
    q12_architecture = dplyr::case_when(
      any(
        gene_class == "circular",
        na.rm = TRUE
      ) ~
        "12q ecDNA",
      
      any(
        gene_class == "noncircular",
        na.rm = TRUE
      ) ~
        "Chromosomal 12q amplification",
      
      TRUE ~
        NA_character_
    ),
    
    .groups = "drop"
  )

q12_patient <- gcap_patient |>
  dplyr::select(
    Patient.ID,
    Tumor_type,
    Time_of_FU,
    event,
    GCAP_class
  ) |>
  dplyr::left_join(
    q12_patient_observed,
    by = "Patient.ID"
  ) |>
  dplyr::mutate(
    q12_any_amplification =
      tidyr::replace_na(
        q12_any_amplification,
        0L
      ),
    
    q12_any_ecDNA =
      tidyr::replace_na(
        q12_any_ecDNA,
        0L
      ),
    
    q12_any_chromosomal =
      tidyr::replace_na(
        q12_any_chromosomal,
        0L
      ),
    
    q12_genes =
      tidyr::replace_na(
        q12_genes,
        ""
      ),
    
    q12_bands =
      tidyr::replace_na(
        q12_bands,
        ""
      ),
    
    q12_architecture =
      tidyr::replace_na(
        q12_architecture,
        "No focal 12q amplification"
      ),
    
    q12_architecture = factor(
      q12_architecture,
      levels = q12_architecture_levels
    ),
    
    Tumor_type = factor(
      Tumor_type,
      levels = c(
        "ARMS",
        "ERMS",
        "SSRMS"
      )
    )
  ) |>
  dplyr::filter(
    !is.na(Tumor_type)
  )


cat(
  "\nPatient-level 12q13-15 architecture:\n"
)

print(
  table(
    q12_patient$Tumor_type,
    q12_patient$q12_architecture,
    useNA = "ifany"
  )
)

# =============================================================================
# 3. BUILD SURVIVAL DATASET
# =============================================================================

q12_survival_data <- q12_patient |>
  dplyr::filter(
    !is.na(Time_of_FU),
    !is.na(event),
    !is.na(q12_architecture)
  )

# =============================================================================
# 4. PAIRWISE FIRTH COX FUNCTION
# =============================================================================

run_q12_firth_contrast <- function(
    data,
    reference_group,
    comparison_group,
    population_label,
    adjust_for_subtype = FALSE
) {
  
  df <- data |>
    dplyr::filter(
      q12_architecture %in%
        c(
          reference_group,
          comparison_group
        )
    ) |>
    dplyr::mutate(
      q12_pair = factor(
        as.character(
          q12_architecture
        ),
        levels = c(
          reference_group,
          comparison_group
        )
      ),
      
      Tumor_type =
        droplevels(
          factor(
            Tumor_type
          )
        )
    ) |>
    droplevels()
  
  reference_n <- sum(
    df$q12_pair == reference_group
  )
  
  comparison_n <- sum(
    df$q12_pair == comparison_group
  )
  
  reference_events <- sum(
    df$event[
      df$q12_pair == reference_group
    ],
    na.rm = TRUE
  )
  
  comparison_events <- sum(
    df$event[
      df$q12_pair == comparison_group
    ],
    na.rm = TRUE
  )
  
  if (
    reference_n < 2 ||
    comparison_n < 2
  ) {
    
    return(
      tibble::tibble(
        Population = population_label,
        Reference = reference_group,
        Comparison = comparison_group,
        Reference_n = reference_n,
        Reference_events = reference_events,
        Comparison_n = comparison_n,
        Comparison_events = comparison_events,
        HR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_,
        Fit_status =
          "Descriptive only: <2 patients in one comparison group"
      )
    )
  }
  
  fit <- tryCatch(
    {
      if (
        adjust_for_subtype &&
        dplyr::n_distinct(
          df$Tumor_type
        ) > 1
      ) {
        
        coxphf::coxphf(
          survival::Surv(
            Time_of_FU,
            event
          ) ~
            q12_pair +
            Tumor_type,
          data = df,
          maxit = 1000,
          maxstep = 0.5
        )
        
      } else {
        
        coxphf::coxphf(
          survival::Surv(
            Time_of_FU,
            event
          ) ~
            q12_pair,
          data = df,
          maxit = 1000,
          maxstep = 0.5
        )
      }
    },
    error = function(e) {
      message(
        "12q Firth Cox failed for ",
        population_label,
        " / ",
        comparison_group,
        " vs ",
        reference_group,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(fit)) {
    
    return(
      tibble::tibble(
        Population = population_label,
        Reference = reference_group,
        Comparison = comparison_group,
        Reference_n = reference_n,
        Reference_events = reference_events,
        Comparison_n = comparison_n,
        Comparison_events = comparison_events,
        HR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_,
        Fit_status = "Firth Cox model failed"
      )
    )
  }
  
  coef_names <- names(
    fit$coefficients
  )
  
  coef_index <- which(
    grepl(
      "^q12_pair",
      coef_names
    )
  )
  
  if (
    length(
      coef_index
    ) != 1
  ) {
    
    return(
      tibble::tibble(
        Population = population_label,
        Reference = reference_group,
        Comparison = comparison_group,
        Reference_n = reference_n,
        Reference_events = reference_events,
        Comparison_n = comparison_n,
        Comparison_events = comparison_events,
        HR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_,
        Fit_status =
          "12q architecture coefficient not found"
      )
    )
  }
  
  i <- coef_index[1]
  
  tibble::tibble(
    Population = population_label,
    Reference = reference_group,
    Comparison = comparison_group,
    Reference_n = reference_n,
    Reference_events = reference_events,
    Comparison_n = comparison_n,
    Comparison_events = comparison_events,
    
    HR = exp(
      fit$coefficients[i]
    ),
    
    # coxphf CI values are already on the HR scale.
    CI_lower =
      fit$ci.lower[i],
    
    CI_upper =
      fit$ci.upper[i],
    
    P_value =
      fit$prob[i],
    
    Fit_status =
      "OK"
  )
}

# =============================================================================
# 5. RUN THREE PRESPECIFIED 12q CONTRASTS
# =============================================================================

run_q12_survival_set <- function(
    data,
    population_label,
    adjust_for_subtype = FALSE
) {
  
  contrasts <- list(
    c(
      "No focal 12q amplification",
      "Chromosomal 12q amplification"
    ),
    c(
      "No focal 12q amplification",
      "12q ecDNA"
    ),
    c(
      "Chromosomal 12q amplification",
      "12q ecDNA"
    )
  )
  
  result <- dplyr::bind_rows(
    lapply(
      contrasts,
      function(x) {
        
        run_q12_firth_contrast(
          data = data,
          reference_group = x[1],
          comparison_group = x[2],
          population_label = population_label,
          adjust_for_subtype =
            adjust_for_subtype
        )
      }
    )
  )
  
  result |>
    dplyr::mutate(
      FDR = stats::p.adjust(
        P_value,
        method = "BH"
      )
    )
}

# Overall analysis adjusted for RMS subtype.
q12_survival_overall_adjusted <-
  run_q12_survival_set(
    data =
      q12_survival_data,
    population_label =
      "Overall — adjusted for RMS subtype",
    adjust_for_subtype =
      TRUE
  )

# Subtype-specific univariable analyses.
q12_survival_by_subtype <-
  dplyr::bind_rows(
    lapply(
      c(
        "ARMS",
        "ERMS",
        "SSRMS"
      ),
      function(t) {
        
        run_q12_survival_set(
          data =
            q12_survival_data |>
            dplyr::filter(
              Tumor_type == t
            ),
          population_label =
            t,
          adjust_for_subtype =
            FALSE
        )
      }
    )
  )

q12_survival_firth_results <-
  dplyr::bind_rows(
    q12_survival_overall_adjusted,
    q12_survival_by_subtype
  )

cat(
  "\n12q architecture-specific Firth Cox survival results:\n"
)

print(
  q12_survival_firth_results
)



# -----------------------------------------------------------------------------
# 19F. ERMS multivariable 12q model
# Reuse the complete clinical variables created in Section 15.
# -----------------------------------------------------------------------------

erms_q12_mv <- erms_multivariable_3group_all |>
  dplyr::select(
    Patient.ID,
    Time_of_FU,
    event,
    TP53_status,
    Age_group,
    Sex,
    Stage_group,
    Clinical_group_binary
  ) |>
  dplyr::inner_join(
    q12_patient |>
      dplyr::filter(
        Tumor_type == "ERMS"
      ) |>
      dplyr::select(
        Patient.ID,
        q12_architecture
      ),
    by = "Patient.ID"
  ) |>
  dplyr::filter(
    !is.na(Time_of_FU),
    !is.na(event),
    !is.na(q12_architecture),
    !is.na(TP53_status),
    !is.na(Age_group),
    !is.na(Sex),
    !is.na(Stage_group),
    !is.na(Clinical_group_binary)
  ) |>
  dplyr::mutate(
    q12_architecture = factor(
      q12_architecture,
      levels = c(
        "No focal 12q amplification",
        "Chromosomal 12q amplification",
        "12q ecDNA"
      )
    )
  ) |>
  droplevels()

cat(
  "\nERMS 12q multivariable cohort: n=",
  nrow(erms_q12_mv),
  "; deaths=",
  sum(erms_q12_mv$event, na.rm = TRUE),
  "\n",
  sep = ""
)

print(
  table(
    erms_q12_mv$q12_architecture,
    erms_q12_mv$event,
    useNA = "ifany"
  )
)

# No focal 12q reference:
#   Chromosomal vs no focal
#   ecDNA vs no focal
fit_erms_q12_mv_no_focal <- coxphf::coxphf(
  survival::Surv(
    Time_of_FU,
    event
  ) ~
    q12_architecture +
    TP53_status +
    Age_group +
    Sex +
    Stage_group +
    Clinical_group_binary,
  data = erms_q12_mv,
  maxit = 1000,
  maxstep = 0.5
)

# Chromosomal 12q reference:
#   ecDNA vs chromosomal
erms_q12_mv_chrom_ref <- erms_q12_mv |>
  dplyr::mutate(
    q12_architecture = stats::relevel(
      q12_architecture,
      ref = "Chromosomal 12q amplification"
    )
  )

fit_erms_q12_mv_chrom <- coxphf::coxphf(
  survival::Surv(
    Time_of_FU,
    event
  ) ~
    q12_architecture +
    TP53_status +
    Age_group +
    Sex +
    Stage_group +
    Clinical_group_binary,
  data = erms_q12_mv_chrom_ref,
  maxit = 1000,
  maxstep = 0.5
)

extract_q12_mv_term <- function(
    fit,
    term_name,
    comparison_label
) {
  
  idx <- match(
    term_name,
    names(fit$coefficients)
  )
  
  if (is.na(idx)) {
    return(
      tibble::tibble(
        Comparison = comparison_label,
        HR = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_
      )
    )
  }
  
  tibble::tibble(
    Comparison = comparison_label,
    HR = exp(fit$coefficients[idx]),
    CI_lower = fit$ci.lower[idx],
    CI_upper = fit$ci.upper[idx],
    P_value = fit$prob[idx]
  )
}

erms_q12_adjusted_results <- dplyr::bind_rows(
  extract_q12_mv_term(
    fit_erms_q12_mv_no_focal,
    "q12_architectureChromosomal 12q amplification",
    "Chromosomal vs no focal"
  ),
  extract_q12_mv_term(
    fit_erms_q12_mv_no_focal,
    "q12_architecture12q ecDNA",
    "ecDNA vs no focal"
  ),
  extract_q12_mv_term(
    fit_erms_q12_mv_chrom,
    "q12_architecture12q ecDNA",
    "ecDNA vs chromosomal"
  )
) |>
  dplyr::mutate(
    FDR = stats::p.adjust(
      P_value,
      method = "BH"
    )
  )

print(erms_q12_adjusted_results)


# -----------------------------------------------------------------------------
# 19G. ARMS / ERMS 12q forest
# ARMS: univariable
# ERMS: univariable + full multivariable
# -----------------------------------------------------------------------------

# 1. Formatting helpers
# -----------------------------------------------------------------------------

fmt_num <- function(x) {
  ifelse(
    is.na(x),
    "NE",
    sub(
      "\\.?0+$",
      "",
      sprintf("%.2f", x)
    )
  )
}

fmt_p <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "NE",
    x < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", x)
  )
}

make_contrast <- function(
    Reference,
    Comparison
) {
  
  dplyr::case_when(
    Reference == "No focal 12q amplification" &
      Comparison == "Chromosomal 12q amplification" ~
      "Chromosomal vs no focal",
    
    Reference == "No focal 12q amplification" &
      Comparison == "12q ecDNA" ~
      "ecDNA vs no focal",
    
    Reference == "Chromosomal 12q amplification" &
      Comparison == "12q ecDNA" ~
      "ecDNA vs chromosomal",
    
    TRUE ~ paste0(
      Comparison,
      " vs ",
      Reference
    )
  )
}

# -----------------------------------------------------------------------------
# 2. ARMS univariable results
# -----------------------------------------------------------------------------

arms_uni <- q12_survival_firth_results |>
  dplyr::filter(
    Population == "ARMS"
  ) |>
  dplyr::mutate(
    Population = "ARMS",
    Contrast = make_contrast(
      Reference,
      Comparison
    ),
    Model = "Univariable"
  ) |>
  dplyr::select(
    Population,
    Contrast,
    Model,
    HR,
    CI_lower,
    CI_upper,
    P_value,
    FDR
  )

# -----------------------------------------------------------------------------
# 3. ERMS univariable results
# -----------------------------------------------------------------------------

erms_uni <- q12_survival_firth_results |>
  dplyr::filter(
    Population == "ERMS"
  ) |>
  dplyr::mutate(
    Population = "ERMS",
    Contrast = make_contrast(
      Reference,
      Comparison
    ),
    Model = "Univariable"
  ) |>
  dplyr::select(
    Population,
    Contrast,
    Model,
    HR,
    CI_lower,
    CI_upper,
    P_value,
    FDR
  )

# -----------------------------------------------------------------------------
# 4. ERMS multivariable results
# -----------------------------------------------------------------------------

erms_mv <- erms_q12_adjusted_results |>
  dplyr::transmute(
    Population = "ERMS",
    Contrast = Comparison,
    Model = "Multivariable",
    HR,
    CI_lower,
    CI_upper,
    P_value,
    FDR
  )

# -----------------------------------------------------------------------------
# 5. Combine
# -----------------------------------------------------------------------------

plot_dat <- dplyr::bind_rows(
  arms_uni,
  erms_uni,
  erms_mv
) |>
  dplyr::mutate(
    HR_CI = dplyr::if_else(
      !is.na(HR),
      paste0(
        fmt_num(HR),
        " (",
        fmt_num(CI_lower),
        "–",
        fmt_num(CI_upper),
        ")"
      ),
      "NE"
    ),
    
    P_text = fmt_p(
      P_value
    ),
    
    Q_text = fmt_p(
      FDR
    ),
    
    significant =
      !is.na(FDR) &
      FDR < 0.05,
    
    estimable =
      !is.na(HR) &
      !is.na(CI_lower) &
      !is.na(CI_upper) &
      HR > 0 &
      CI_lower > 0 &
      CI_upper > 0
  )

# -----------------------------------------------------------------------------
# 6. Fixed row positions
#
# ARMS:
#   3 comparisons × 1 row
#
# ERMS:
#   3 comparisons × 2 rows
# -----------------------------------------------------------------------------

row_map <- data.frame(
  Population = c(
    "ARMS",
    "ARMS",
    "ARMS",
    
    "ERMS",
    "ERMS",
    "ERMS",
    "ERMS",
    "ERMS",
    "ERMS"
  ),
  
  Contrast = c(
    "Chromosomal vs no focal",
    "ecDNA vs no focal",
    "ecDNA vs chromosomal",
    
    "Chromosomal vs no focal",
    "Chromosomal vs no focal",
    
    "ecDNA vs no focal",
    "ecDNA vs no focal",
    
    "ecDNA vs chromosomal",
    "ecDNA vs chromosomal"
  ),
  
  Model = c(
    "Univariable",
    "Univariable",
    "Univariable",
    
    "Univariable",
    "Multivariable",
    
    "Univariable",
    "Multivariable",
    
    "Univariable",
    "Multivariable"
  ),
  
  y = c(
    12.0,
    10.8,
    9.6,
    
    6.8,
    6.0,
    
    4.4,
    3.6,
    
    2.0,
    1.2
  ),
  
  stringsAsFactors = FALSE
)

plot_dat <- plot_dat |>
  dplyr::left_join(
    row_map,
    by = c(
      "Population",
      "Contrast",
      "Model"
    )
  ) |>
  dplyr::arrange(
    dplyr::desc(y)
  )

# -----------------------------------------------------------------------------
# 7. Section and comparison headers
# -----------------------------------------------------------------------------

population_headers <- data.frame(
  Population = c(
    "ARMS",
    "ERMS"
  ),
  
  y = c(
    13.15,
    8.15
  ),
  
  stringsAsFactors = FALSE
)

contrast_headers <- data.frame(
  Population = c(
    "ARMS",
    "ARMS",
    "ARMS",
    "ERMS",
    "ERMS",
    "ERMS"
  ),
  
  Contrast = c(
    "Chromosomal vs no focal",
    "ecDNA vs no focal",
    "ecDNA vs chromosomal",
    "Chromosomal vs no focal",
    "ecDNA vs no focal",
    "ecDNA vs chromosomal"
  ),
  
  y = c(
    12.55,
    11.35,
    10.15,
    7.35,
    4.95,
    2.55
  ),
  
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 8. Gray comparison bands
#
# ARMS:
#   gray / white / gray
#
# ERMS:
#   gray / white / gray
# -----------------------------------------------------------------------------

band_dat <- data.frame(
  ymin = c(
    11.75,
    9.35,
    5.65,
    0.85
  ),
  
  ymax = c(
    12.80,
    10.40,
    7.65,
    2.75
  )
)

# -----------------------------------------------------------------------------

# =============================================================================
# IMPROVED FOREST PLOT
#
# Main changes:
#   1. Exact shared y-axis limits between table and forest
#   2. cowplot::align_plots() forces plotting panels to line up vertically
#   3. No significance-based bolding
#   4. Same point size / line width for every estimate
#   5. Cleaner spacing between table and forest
# =============================================================================

library(cowplot)


# -----------------------------------------------------------------------------
# Common vertical plotting range
# -----------------------------------------------------------------------------

y_limits <- c(
  0.70,
  14.55
)


# -----------------------------------------------------------------------------
# 9. LEFT TABLE
# -----------------------------------------------------------------------------

p_table <- ggplot2::ggplot() +
  
  # Alternating comparison bands
  ggplot2::geom_rect(
    data = band_dat,
    ggplot2::aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "grey94",
    color = NA,
    inherit.aes = FALSE
  ) +
  
  # ARMS / ERMS headers
  ggplot2::geom_text(
    data = population_headers,
    ggplot2::aes(
      x = 0,
      y = y,
      label = Population
    ),
    hjust = 0,
    fontface = "bold",
    size = 5.1
  ) +
  
  # Comparison headers
  ggplot2::geom_text(
    data = contrast_headers,
    ggplot2::aes(
      x = 0.02,
      y = y,
      label = Contrast
    ),
    hjust = 0,
    fontface = "bold",
    size = 4.25
  ) +
  
  # Model
  ggplot2::geom_text(
    data = plot_dat,
    ggplot2::aes(
      x = 0.06,
      y = y,
      label = Model
    ),
    hjust = 0,
    size = 3.9
  ) +
  
  # HR and 95% CI
  ggplot2::geom_text(
    data = plot_dat,
    ggplot2::aes(
      x = 0.51,
      y = y,
      label = HR_CI
    ),
    hjust = 0.5,
    size = 3.9
  ) +
  
  # P value
  ggplot2::geom_text(
    data = plot_dat,
    ggplot2::aes(
      x = 0.80,
      y = y,
      label = P_text
    ),
    hjust = 0.5,
    size = 3.9
  ) +
  
  # BH q
  ggplot2::geom_text(
    data = plot_dat,
    ggplot2::aes(
      x = 0.95,
      y = y,
      label = Q_text
    ),
    hjust = 0.5,
    size = 3.9
  ) +
  
  # Column headers
  ggplot2::annotate(
    "text",
    x = c(
      0.06,
      0.51,
      0.80,
      0.95
    ),
    y = 14.15,
    label = c(
      "Model",
      "HR (95% CI)",
      "P",
      "BH q"
    ),
    hjust = c(
      0,
      0.5,
      0.5,
      0.5
    ),
    fontface = "bold",
    size = 4.25
  ) +
  
  # IMPORTANT:
  # Exact same y scale as forest panel
  ggplot2::scale_y_continuous(
    limits = y_limits,
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      1.03
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0)
    )
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 5,
      r = 8,
      b = 5,
      l = 8
    )
  )


# -----------------------------------------------------------------------------
# 10. RIGHT FOREST
# -----------------------------------------------------------------------------

forest_dat <- plot_dat |>
  dplyr::filter(
    estimable
  )


p_forest <- ggplot2::ggplot() +
  
  # Same alternating comparison bands
  ggplot2::geom_rect(
    data = band_dat,
    ggplot2::aes(
      xmin = 0.1,
      xmax = 100,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "grey94",
    color = NA,
    inherit.aes = FALSE
  ) +
  
  # Null HR
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  # Confidence intervals
  ggplot2::geom_segment(
    data = forest_dat,
    ggplot2::aes(
      x = CI_lower,
      xend = CI_upper,
      y = y,
      yend = y
    ),
    linewidth = 0.8,
    lineend = "butt"
  ) +
  
  # HR estimates
  ggplot2::geom_point(
    data = forest_dat,
    ggplot2::aes(
      x = HR,
      y = y,
      shape = Model
    ),
    size = 3.4
  ) +
  
  # Non-estimable rows
  ggplot2::geom_text(
    data = plot_dat |>
      dplyr::filter(
        !estimable
      ),
    ggplot2::aes(
      x = 1,
      y = y,
      label = "NE"
    ),
    size = 3.8
  ) +
  
  # Square = univariable
  # Diamond = multivariable
  ggplot2::scale_shape_manual(
    values = c(
      "Univariable" = 15,
      "Multivariable" = 18
    )
  ) +
  
  # Log HR axis
  ggplot2::scale_x_log10(
    limits = c(
      0.1,
      100
    ),
    breaks = c(
      0.1,
      0.5,
      1,
      5,
      20,
      100
    ),
    labels = c(
      "0.1",
      "0.5",
      "1",
      "5",
      "20",
      "100"
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  # IMPORTANT:
  # Identical y scale to table
  ggplot2::scale_y_continuous(
    limits = y_limits,
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  
  ggplot2::theme_classic(
    base_size = 12
  ) +
  
  ggplot2::theme(
    
    axis.text.y =
      ggplot2::element_blank(),
    
    axis.ticks.y =
      ggplot2::element_blank(),
    
    axis.title.y =
      ggplot2::element_blank(),
    
    axis.line.y =
      ggplot2::element_blank(),
    
    axis.text.x =
      ggplot2::element_text(
        size = 10.5
      ),
    
    axis.title.x =
      ggplot2::element_text(
        size = 11.5,
        margin = ggplot2::margin(
          t = 7
        )
      ),
    
    legend.position = "none",
    
    plot.margin =
      ggplot2::margin(
        t = 5,
        r = 8,
        b = 5,
        l = 4
      )
  ) +
  
  ggplot2::labs(
    x = "Hazard ratio"
  )


# -----------------------------------------------------------------------------
# 11. FORCE TABLE AND FOREST PANELS TO BE VERTICALLY ALIGNED
# -----------------------------------------------------------------------------

aligned_panels <- cowplot::align_plots(
  p_table,
  p_forest,
  align = "v",
  axis = "tb"
)


forest_body <- cowplot::plot_grid(
  aligned_panels[[1]],
  aligned_panels[[2]],
  nrow = 1,
  rel_widths = c(
    3.15,
    1.55
  )
)


# -----------------------------------------------------------------------------
# Title
# -----------------------------------------------------------------------------

forest_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "12q13–15 amplification architecture and overall survival",
    x = 0.5,
    hjust = 0.5,
    fontface = "bold",
    size = 16
  )


# -----------------------------------------------------------------------------
# Footnote
# -----------------------------------------------------------------------------

forest_note <- cowplot::ggdraw() +
  cowplot::draw_label(
    paste0(
      "Subtype-specific Firth Cox regression. ",
      "ARMS estimates are univariable. ",
      "ERMS multivariable models adjust for TP53 status, age, sex, ",
      "local stage, and clinical group. ",
      "BH q values were calculated across the three architecture contrasts ",
      "within each analysis."
    ),
    x = 0.01,
    hjust = 0,
    vjust = 0.5,
    size = 8.5
  )


# -----------------------------------------------------------------------------
# Final figure
# -----------------------------------------------------------------------------

arms_erms_q12_forest <- cowplot::plot_grid(
  forest_title,
  forest_body,
  forest_note,
  ncol = 1,
  rel_heights = c(
    0.08,
    1,
    0.08
  )
)


# Display
print(
  arms_erms_q12_forest
)


# -----------------------------------------------------------------------------
# Display forest plot
# -----------------------------------------------------------------------------

