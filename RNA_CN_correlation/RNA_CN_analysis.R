library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(readr)
library(lme4)

RNA <- file.path("data", "PURPLE_geneCN_Salmon_TPM_matched.tsv")
amplicon_file <- file.path("data", "Joined_final.csv")


genes_to_test <- c(
  "MYC",
  "MYCN",
  "MDM2",
  "CDK4",
  "GLI1", 
  "FOXO1", 
  "PAX7", 
  "DDIT3", 
  "FGFR1", "YEATS4"
)


architecture_colors <- c(
  "NoFocalAmp" = "lightgrey",
  "ChrAmp"      = "#EE0000FF",
  "ecDNA"       = "#6C77AD"
)

architecture_labels <- c(
  "NoFocalAmp" = "No focal amp",
  "ChrAmp"      = "Chromosomal",
  "ecDNA"       = "ecDNA"
)


# Detailed architecture colors used ONLY for the continuous CN-RNA plot.
continuous_architecture_colors <- c(
  "NoFocalAmp" = "lightgrey",
  "Linear" = "#F0E685",
  "Complex-non-cyclic" = "#CE3D32",
  "BFB" = "#CC79A7",
  "ecDNA" = "#6C77AD"
)

continuous_architecture_labels <- c(
  "NoFocalAmp" = "No focal amp",
  "Linear" = "Linear",
  "Complex-non-cyclic" = "Complex non-cyclic",
  "BFB" = "BFB",
  "ecDNA" = "ecDNA"
)



# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

extract_genes <- function(x) {
  
  if (
    is.na(x) ||
    x == "" ||
    x == "[]" ||
    x == "NA"
  ) {
    return(character(0))
  }
  
  x <- gsub(
    "\\[|\\]|'",
    "",
    x
  )
  
  genes <- trimws(
    unlist(
      strsplit(
        x,
        ","
      )
    )
  )
  
  genes[genes != ""]
}


clean_amplicon_class <- function(x) {
  
  dplyr::case_when(
    
    stringr::str_detect(
      x,
      stringr::regex(
        "^ecDNA$",
        ignore_case = TRUE
      )
    ) ~ "ecDNA",
    
    stringr::str_detect(
      x,
      stringr::regex(
        "BFB",
        ignore_case = TRUE
      )
    ) ~ "BFB",
    
    stringr::str_detect(
      x,
      stringr::regex(
        "Linear",
        ignore_case = TRUE
      )
    ) ~ "Linear",
    
    stringr::str_detect(
      x,
      stringr::regex(
        "Complex",
        ignore_case = TRUE
      )
    ) ~ "Complex",
    
    TRUE ~ "Other"
  )
}


spearman_by_gene <- function(
    data,
    cn_col,
    amplified_only = FALSE
) {
  
  cn_col <- rlang::ensym(cn_col)
  
  x <- data |>
    dplyr::filter(
      !is.na(!!cn_col),
      !is.na(TPM)
    )
  
  if (amplified_only) {
    
    x <- x |>
      dplyr::filter(
        amplified
      )
  }
  
  x |>
    dplyr::group_by(
      gene
    ) |>
    dplyr::group_modify(
      \(.x, .y) {
        
        cn <- dplyr::pull(
          .x,
          !!cn_col
        )
        
        if (
          nrow(.x) < 3 ||
          dplyr::n_distinct(cn) < 2
        ) {
          
          return(
            tibble::tibble(
              n = nrow(.x),
              rho = NA_real_,
              p_value = NA_real_
            )
          )
        }
        
        test <- suppressWarnings(
          stats::cor.test(
            cn,
            .x$TPM,
            method = "spearman",
            exact = FALSE
          )
        )
        
        tibble::tibble(
          n = nrow(.x),
          rho = unname(
            test$estimate
          ),
          p_value = test$p.value
        )
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      FDR = stats::p.adjust(
        p_value,
        method = "BH"
      )
    ) |>
    dplyr::arrange(
      FDR
    )
}

# =============================================================================
# 3. READ MATCHED PURPLE + SALMON DATA
# =============================================================================

cn_rna <- read.delim(
  RNA,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cn_rna <- c(
  "sample_name",
  "gene",
  "minCopyNumber",
  "maxCopyNumber",
  "TPM"
)

stopifnot(
  all(
    required_cn_rna %in%
      names(cn_rna)
  )
)

cn_rna <- cn_rna |>
  dplyr::mutate(
    minCopyNumber = as.numeric(
      minCopyNumber
    ),
    maxCopyNumber = as.numeric(
      maxCopyNumber
    ),
    TPM = as.numeric(
      TPM
    ),
    log2_minCN = log2(
      minCopyNumber + 1
    ),
    log2_maxCN = log2(
      maxCopyNumber + 1
    ),
    log2_TPM = log2(
      TPM + 1
    ),
    CN_span =
      maxCopyNumber -
      minCopyNumber
  ) |>
  dplyr::filter(
    gene %in%
      genes_to_test
  ) |>
  dplyr::distinct(
    sample_name,
    gene,
    .keep_all = TRUE
  )

# =============================================================================
# 4. READ AMPLICONSUITE TABLE
# =============================================================================

amp <- read.csv(
  amplicon_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Patient identifier is required for patient-aware sensitivity analyses.
if (!"patient" %in% names(amp)) {
  stop("The AmpliconSuite table must contain a 'patient' column.")
}

sample_patient_map <- amp |>
  dplyr::filter(!is.na(Sample.name), !is.na(patient), patient != "") |>
  dplyr::distinct(sample_name = Sample.name, patient)

if (sample_patient_map |> dplyr::count(sample_name) |> dplyr::filter(n > 1) |> nrow() > 0) {
  stop("At least one sample maps to more than one patient.")
}

bad_cols <- which(
  is.na(names(amp)) |
    names(amp) == ""
)

if (
  length(bad_cols) > 0
) {
  
  amp <- amp[
    ,
    -bad_cols,
    drop = FALSE
  ]
}

names(amp) <- make.unique(
  names(amp)
)

required_amp <- c(
  "Sample.name",
  "Classification",
  "All.genes"
)

stopifnot(
  all(
    required_amp %in%
      names(amp)
  )
)

amp <- amp |>
  dplyr::mutate(
    Classification_clean =
      clean_amplicon_class(
        Classification
      )
  )

# =============================================================================
# 5. CONVERT AMPLICON FEATURES TO SAMPLE-GENE ROWS
# =============================================================================

amp_gene_long <- amp |>
  dplyr::mutate(
    gene_list = purrr::map(
      `All.genes`,
      extract_genes
    )
  ) |>
  dplyr::select(
    sample_name = Sample.name,
    Classification =
      Classification_clean,
    gene_list
  ) |>
  tidyr::unnest_longer(
    gene_list,
    values_to = "gene"
  ) |>
  dplyr::filter(
    gene %in%
      genes_to_test,
    Classification %in%
      c(
        "ecDNA",
        "BFB",
        "Linear",
        "Complex"
      )
  )

# =============================================================================
# 6. ASSIGN GENE-LEVEL AMPLIFICATION ARCHITECTURE
# =============================================================================

gene_architecture <- amp_gene_long |>
  dplyr::group_by(
    sample_name,
    gene
  ) |>
  dplyr::summarise(
    
    has_ecDNA = any(
      Classification ==
        "ecDNA"
    ),
    
    has_ChrAmp = any(
      Classification %in%
        c(
          "BFB",
          "Linear",
          "Complex"
        )
    ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    architecture =
      dplyr::case_when(
        
        has_ecDNA &
          has_ChrAmp ~
          "Mixed",
        
        has_ecDNA ~
          "ecDNA",
        
        has_ChrAmp ~
          "ChrAmp",
        
        TRUE ~
          NA_character_
      )
  ) |>
  dplyr::select(
    sample_name,
    gene,
    architecture
  )


# Detailed architecture used ONLY for the continuous CN-RNA scatter plot.
# Boxplots and all statistical analyses continue to use the broad
# NoFocalAmp / ChrAmp / ecDNA classification above.
gene_architecture_detailed <- amp_gene_long |>
  dplyr::group_by(
    sample_name,
    gene
  ) |>
  dplyr::summarise(
    has_ecDNA = any(Classification == "ecDNA"),
    has_BFB = any(Classification == "BFB"),
    has_Complex = any(Classification == "Complex"),
    has_Linear = any(Classification == "Linear"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    architecture_detailed = dplyr::case_when(
      has_ecDNA ~ "ecDNA",
      has_BFB ~ "BFB",
      has_Complex ~ "Complex-non-cyclic",
      has_Linear ~ "Linear",
      TRUE ~ NA_character_
    )
  ) |>
  dplyr::select(
    sample_name,
    gene,
    architecture_detailed
  )

continuous_plot_data <- cn_rna |>
  dplyr::left_join(
    gene_architecture_detailed,
    by = c("sample_name", "gene")
  ) |>
  dplyr::mutate(
    architecture_detailed = tidyr::replace_na(
      architecture_detailed,
      "NoFocalAmp"
    ),
    architecture_detailed = factor(
      architecture_detailed,
      levels = c(
        "NoFocalAmp",
        "Linear",
        "Complex-non-cyclic",
        "BFB",
        "ecDNA"
      )
    )
  )

# =============================================================================
# 7. CREATE FINAL ANALYSIS DATASET
# =============================================================================

analysis_data <- cn_rna |>
  dplyr::left_join(
    gene_architecture,
    by = c(
      "sample_name",
      "gene"
    )
  ) |>
  dplyr::left_join(
    sample_patient_map,
    by = "sample_name"
  ) |>
  dplyr::mutate(
    
    architecture =
      tidyr::replace_na(
        architecture,
        "NoFocalAmp"
      )
  ) |>
  dplyr::filter(
    architecture %in%
      c(
        "NoFocalAmp",
        "ChrAmp",
        "ecDNA"
      )
  ) |>
  dplyr::mutate(
    
    architecture = factor(
      architecture,
      levels = c(
        "NoFocalAmp",
        "ChrAmp",
        "ecDNA"
      )
    ),
    
    amplified =
      architecture !=
      "NoFocalAmp"
  )

# =============================================================================
# 8. GROUP COUNTS
# =============================================================================

group_counts <- analysis_data |>
  dplyr::group_by(
    gene,
    architecture
  ) |>
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  )

# =============================================================================
# 9. AMPLIFIED VS NON-AMPLIFIED RNA EXPRESSION
# =============================================================================

amp_vs_noamp <- analysis_data |>
  dplyr::filter(
    !is.na(log2_TPM)
  ) |>
  dplyr::group_by(
    gene
  ) |>
  dplyr::group_modify(
    \(.x, .y) {
      
      n_amp <- sum(
        .x$amplified
      )
      
      n_noamp <- sum(
        !.x$amplified
      )
      
      p_value <- NA_real_
      
      if (
        n_amp >= 2 &&
        n_noamp >= 2
      ) {
        
        p_value <-
          stats::wilcox.test(
            log2_TPM ~ amplified,
            data = .x,
            exact = FALSE
          )$p.value
      }
      
      tibble::tibble(
        
        n_total =
          nrow(.x),
        
        n_amplified =
          n_amp,
        
        n_non_amplified =
          n_noamp,
        
        median_log2TPM_amplified =
          if (n_amp > 0) {
            median(
              .x$log2_TPM[
                .x$amplified
              ],
              na.rm = TRUE
            )
          } else {
            NA_real_
          },
        
        median_log2TPM_non_amplified =
          if (n_noamp > 0) {
            median(
              .x$log2_TPM[
                !.x$amplified
              ],
              na.rm = TRUE
            )
          } else {
            NA_real_
          },
        
        p_value =
          p_value
      )
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR = stats::p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  dplyr::arrange(
    FDR
  )

# =============================================================================
# 10. THREE-GROUP ARCHITECTURE TEST
#
# NoFocalAmp vs ChrAmp vs ecDNA
# Overall Kruskal-Wallis test.
# =============================================================================

architecture_test <- analysis_data |>
  dplyr::filter(
    !is.na(log2_TPM)
  ) |>
  dplyr::group_by(
    gene
  ) |>
  dplyr::group_modify(
    \(.x, .y) {
      
      n_groups <-
        dplyr::n_distinct(
          .x$architecture[
            !is.na(
              .x$architecture
            )
          ]
        )
      
      p_value <- NA_real_
      
      if (
        n_groups >= 2 &&
        nrow(.x) >= 3
      ) {
        
        p_value <-
          stats::kruskal.test(
            log2_TPM ~ architecture,
            data = .x
          )$p.value
      }
      
      tibble::tibble(
        p_value =
          p_value
      )
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    
    FDR =
      stats::p.adjust(
        p_value,
        method = "BH"
      ),
    
    FDR_label =
      dplyr::case_when(
        
        is.na(FDR) ~
          "",
        
        FDR < 0.001 ~
          "FDR < 0.001",
        
        TRUE ~
          paste0(
            "FDR = ",
            formatC(
              FDR,
              format = "f",
              digits = 3
            )
          )
      )
  )

# =============================================================================
# 11. CONTINUOUS PURPLE MINIMUM CN VS RNA
# =============================================================================

cn_rna_all <- spearman_by_gene(
  analysis_data,
  minCopyNumber
)

# =============================================================================
# 12. CONTINUOUS CN VS RNA AMONG AMPLIFIED CASES
# =============================================================================

cn_rna_amplified <- spearman_by_gene(
  analysis_data,
  minCopyNumber,
  amplified_only = TRUE
) |>
  dplyr::transmute(
    gene,
    n_amplified = n,
    rho_amplified = rho,
    p_amplified = p_value,
    FDR_amplified = FDR
  )

# =============================================================================
# 13. MAXIMUM COPY NUMBER SENSITIVITY ANALYSIS
# =============================================================================

cn_rna_max <- spearman_by_gene(
  analysis_data,
  maxCopyNumber
) |>
  dplyr::transmute(
    gene,
    n_maxCN = n,
    rho_maxCN = rho,
    p_maxCN = p_value,
    FDR_maxCN = FDR
  )

# =============================================================================
# 14. ecDNA VS CHROMOSOMAL AMPLIFICATION
# =============================================================================

ecdna_vs_chramp <- analysis_data |>
  dplyr::filter(
    architecture %in%
      c(
        "ChrAmp",
        "ecDNA"
      ),
    !is.na(log2_TPM)
  ) |>
  droplevels() |>
  dplyr::group_by(
    gene
  ) |>
  dplyr::group_modify(
    \(.x, .y) {
      
      n_ecDNA <- sum(
        .x$architecture ==
          "ecDNA"
      )
      
      n_ChrAmp <- sum(
        .x$architecture ==
          "ChrAmp"
      )
      
      p_value <- NA_real_
      
      if (
        n_ecDNA >= 2 &&
        n_ChrAmp >= 2
      ) {
        
        p_value <-
          stats::wilcox.test(
            log2_TPM ~ architecture,
            data = .x,
            exact = FALSE
          )$p.value
      }
      
      tibble::tibble(
        
        n_ecDNA =
          n_ecDNA,
        
        n_ChrAmp =
          n_ChrAmp,
        
        median_log2TPM_ecDNA =
          if (n_ecDNA > 0) {
            median(
              .x$log2_TPM[
                .x$architecture ==
                  "ecDNA"
              ],
              na.rm = TRUE
            )
          } else {
            NA_real_
          },
        
        median_log2TPM_ChrAmp =
          if (n_ChrAmp > 0) {
            median(
              .x$log2_TPM[
                .x$architecture ==
                  "ChrAmp"
              ],
              na.rm = TRUE
            )
          } else {
            NA_real_
          },
        
        p_value =
          p_value
      )
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    FDR =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  ) |>
  dplyr::arrange(
    FDR
  )

# =============================================================================
# 15. COMBINED SUMMARY TABLE
# =============================================================================

summary_results <- amp_vs_noamp |>
  
  dplyr::left_join(
    
    cn_rna_all |>
      dplyr::transmute(
        gene,
        n_CN_all = n,
        rho_CN_all = rho,
        p_CN_all = p_value,
        FDR_CN_all = FDR
      ),
    
    by = "gene"
  ) |>
  
  dplyr::left_join(
    cn_rna_amplified,
    by = "gene"
  ) |>
  
  dplyr::left_join(
    
    cn_rna_max |>
      dplyr::select(
        gene,
        rho_maxCN,
        p_maxCN,
        FDR_maxCN
      ),
    
    by = "gene"
  ) |>
  
  dplyr::left_join(
    
    ecdna_vs_chramp |>
      dplyr::transmute(
        gene,
        n_ecDNA,
        n_ChrAmp,
        median_log2TPM_ecDNA,
        median_log2TPM_ChrAmp,
        p_ecDNA_vs_ChrAmp =
          p_value,
        FDR_ecDNA_vs_ChrAmp =
          FDR
      ),
    
    by = "gene"
  ) |>
  
  dplyr::left_join(
    
    architecture_test |>
      dplyr::transmute(
        gene,
        p_architecture =
          p_value,
        FDR_architecture =
          FDR
      ),
    
    by = "gene"
  )


# =============================================================================
# 17. RNA expression by amplification architecture.
# FDR = overall Kruskal-Wallis FDR for each gene.
# =============================================================================

figure1_stats <- analysis_data |>
  dplyr::group_by(
    gene
  ) |>
  dplyr::summarise(
    y_position =
      max(
        log2_TPM,
        na.rm = TRUE
      ) + 0.5,
    .groups = "drop"
  ) |>
  dplyr::left_join(
    
    architecture_test |>
      dplyr::select(
        gene,
        FDR_label
      ),
    
    by = "gene"
  )

p_expression <- ggplot2::ggplot(
  
  analysis_data,
  
  ggplot2::aes(
    x = architecture,
    y = log2_TPM,
    fill = architecture
  )
) +
  
  ggplot2::geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  
  ggplot2::geom_jitter(
    ggplot2::aes(
      color = architecture
    ),
    width = 0.12,
    size = 1.8,
    alpha = 0.8
  ) +
  
  ggplot2::geom_text(
    data = figure1_stats,
    ggplot2::aes(
      x = 2,
      y = y_position,
      label = FDR_label
    ),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  
  ggplot2::facet_wrap(
    ~ gene,
    scales = "free_y"
  ) +
  
  ggplot2::scale_fill_manual(
    values =
      architecture_colors,
    labels =
      architecture_labels
  ) +
  
  ggplot2::scale_color_manual(
    values =
      architecture_colors,
    labels =
      architecture_labels
  ) +
  
  ggplot2::scale_x_discrete(
    labels =
      architecture_labels
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = expression(
      log[2] *
        "(TPM + 1)"
    ),
    fill = "Architecture",
    color = "Architecture"
  ) +
  
  ggplot2::theme_classic(
    base_size = 12
  ) +
  
  ggplot2::theme(
    
    legend.position =
      "bottom",
    
    axis.text.x =
      ggplot2::element_text(
        angle = 30,
        hjust = 1
      ),
    
    strip.background =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        face = "bold"
      )
  )

)

# =============================================================================
# 18.PURPLE minimum gene copy number vs RNA expression.
# FDR = Spearman correlation FDR for each gene.
# =============================================================================

cn_plot_stats <- cn_rna_all |>
  dplyr::mutate(
    
    stat_label =
      dplyr::case_when(
        
        is.na(rho) ~
          "",
        
        FDR < 0.001 ~
          paste0(
            "\u03c1 = ",
            sprintf(
              "%.2f",
              rho
            ),
            "\nFDR < 0.001"
          ),
        
        TRUE ~
          paste0(
            "\u03c1 = ",
            sprintf(
              "%.2f",
              rho
            ),
            "\nFDR = ",
            sprintf(
              "%.3f",
              FDR
            )
          )
      )
  )

p_cn <- ggplot2::ggplot(
  continuous_plot_data,
  ggplot2::aes(
    x = log2_minCN,
    y = log2_TPM,
    color = architecture_detailed
  )
) +
  
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.85
  ) +
  
  ggplot2::geom_smooth(
    ggplot2::aes(
      group = 1
    ),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.6
  ) +
  
  ggplot2::geom_text(
    data = cn_plot_stats,
    ggplot2::aes(
      x = -Inf,
      y = Inf,
      label = stat_label
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.1,
    size = 3.4
  ) +
  
  ggplot2::facet_wrap(
    ~ gene,
    scales = "free"
  ) +
  
  ggplot2::scale_color_manual(
    values = continuous_architecture_colors,
    labels = continuous_architecture_labels,
    drop = FALSE
  ) +
  
  ggplot2::labs(
    x = expression(
      log[2] *
        "(PURPLE minimum copy number + 1)"
    ),
    y = expression(
      log[2] *
        "(TPM + 1)"
    ),
    color = "Architecture"
  ) +
  
  ggplot2::theme_classic(
    base_size = 12
  ) +
  
  ggplot2::theme(
    legend.position = "bottom",
    
    strip.background =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        face = "bold"
      )
  )

# =============================================================================
# 19. PATIENT-AWARE SENSITIVITY ANALYSES
# =============================================================================
# -----------------------------------------------------------------------------
# 19A. Repeated patients
# -----------------------------------------------------------------------------

patient_sample_counts <- analysis_data |>
  dplyr::filter(!is.na(patient), patient != "") |>
  dplyr::distinct(patient, sample_name) |>
  dplyr::count(patient, name = "n_samples") |>
  dplyr::arrange(dplyr::desc(n_samples), patient)

repeated_patients <- patient_sample_counts |>
  dplyr::filter(n_samples > 1)


# -----------------------------------------------------------------------------
# 19B. One-sample-per-patient sensitivity analysis
# -----------------------------------------------------------------------------

representative_samples <- analysis_data |>
  dplyr::filter(!is.na(patient), patient != "") |>
  dplyr::distinct(patient, sample_name) |>
  dplyr::arrange(patient, sample_name) |>
  dplyr::group_by(patient) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup()

analysis_data_one_per_patient <- analysis_data |>
  dplyr::inner_join(representative_samples, by = c("patient", "sample_name"))

amp_vs_noamp_one_per_patient <- analysis_data_one_per_patient |>
  dplyr::filter(!is.na(log2_TPM)) |>
  dplyr::group_by(gene) |>
  dplyr::group_modify(\(.x, .y) {
    n_amp <- sum(.x$amplified)
    n_noamp <- sum(!.x$amplified)
    p_value <- NA_real_
    if (n_amp >= 2 && n_noamp >= 2) {
      p_value <- stats::wilcox.test(log2_TPM ~ amplified, data = .x, exact = FALSE)$p.value
    }
    tibble::tibble(n_total = nrow(.x), n_amplified = n_amp, n_non_amplified = n_noamp, p_value = p_value)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(FDR = stats::p.adjust(p_value, method = "BH"))

architecture_test_one_per_patient <- analysis_data_one_per_patient |>
  dplyr::filter(!is.na(log2_TPM)) |>
  dplyr::group_by(gene) |>
  dplyr::group_modify(\(.x, .y) {
    n_groups <- dplyr::n_distinct(.x$architecture[!is.na(.x$architecture)])
    p_value <- NA_real_
    if (n_groups >= 2 && nrow(.x) >= 3) {
      p_value <- stats::kruskal.test(log2_TPM ~ architecture, data = .x)$p.value
    }
    tibble::tibble(n = nrow(.x), p_value = p_value)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(FDR = stats::p.adjust(p_value, method = "BH"))

ecdna_vs_chramp_one_per_patient <- analysis_data_one_per_patient |>
  dplyr::filter(architecture %in% c("ChrAmp", "ecDNA"), !is.na(log2_TPM)) |>
  droplevels() |>
  dplyr::group_by(gene) |>
  dplyr::group_modify(\(.x, .y) {
    n_ecDNA <- sum(.x$architecture == "ecDNA")
    n_ChrAmp <- sum(.x$architecture == "ChrAmp")
    p_value <- NA_real_
    if (n_ecDNA >= 2 && n_ChrAmp >= 2) {
      p_value <- stats::wilcox.test(log2_TPM ~ architecture, data = .x, exact = FALSE)$p.value
    }
    tibble::tibble(n_ecDNA = n_ecDNA, n_ChrAmp = n_ChrAmp, p_value = p_value)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(FDR = stats::p.adjust(p_value, method = "BH"))

cn_rna_one_per_patient <- spearman_by_gene(
  analysis_data_one_per_patient,
  minCopyNumber
)


# -----------------------------------------------------------------------------
# 19C. Mixed-effects architecture model
# log2(TPM + 1) ~ architecture + (1 | patient)
# -----------------------------------------------------------------------------

mixed_architecture_test <- analysis_data |>
  dplyr::filter(!is.na(patient), patient != "", !is.na(log2_TPM)) |>
  dplyr::group_by(gene) |>
  dplyr::group_modify(\(.x, .y) {
    x <- droplevels(.x)
    n_groups <- dplyr::n_distinct(x$architecture)
    n_patients <- dplyr::n_distinct(x$patient)
    if (nrow(x) < 4 || n_groups < 2 || n_patients < 2) {
      return(tibble::tibble(n = nrow(x), n_patients = n_patients, chi_square = NA_real_, df = NA_real_, p_value = NA_real_, singular = NA))
    }
    fit_null <- tryCatch(lme4::lmer(log2_TPM ~ 1 + (1 | patient), data = x, REML = FALSE), error = function(e) NULL)
    fit_full <- tryCatch(lme4::lmer(log2_TPM ~ architecture + (1 | patient), data = x, REML = FALSE), error = function(e) NULL)
    if (is.null(fit_null) || is.null(fit_full)) {
      return(tibble::tibble(n = nrow(x), n_patients = n_patients, chi_square = NA_real_, df = NA_real_, p_value = NA_real_, singular = NA))
    }
    lrt <- stats::anova(fit_null, fit_full)
    tibble::tibble(
      n = nrow(x),
      n_patients = n_patients,
      chi_square = lrt[2, "Chisq"],
      df = lrt[2, "Df"],
      p_value = lrt[2, "Pr(>Chisq)"],
      singular = lme4::isSingular(fit_full, tol = 1e-4)
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(FDR = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(FDR)

# -----------------------------------------------------------------------------
# 19D. Mixed-effects ecDNA vs chromosomal amplification
# -----------------------------------------------------------------------------

mixed_ecdna_vs_chramp <- analysis_data |>
  dplyr::filter(architecture %in% c("ChrAmp", "ecDNA"), !is.na(patient), patient != "", !is.na(log2_TPM)) |>
  droplevels() |>
  dplyr::group_by(gene) |>
  dplyr::group_modify(\(.x, .y) {
    x <- droplevels(.x)
    n_ecDNA <- sum(x$architecture == "ecDNA")
    n_ChrAmp <- sum(x$architecture == "ChrAmp")
    n_patients <- dplyr::n_distinct(x$patient)
    if (n_ecDNA < 2 || n_ChrAmp < 2 || n_patients < 2) {
      return(tibble::tibble(n_ecDNA = n_ecDNA, n_ChrAmp = n_ChrAmp, n_patients = n_patients, chi_square = NA_real_, df = NA_real_, p_value = NA_real_, singular = NA))
    }
    fit_null <- tryCatch(lme4::lmer(log2_TPM ~ 1 + (1 | patient), data = x, REML = FALSE), error = function(e) NULL)
    fit_full <- tryCatch(lme4::lmer(log2_TPM ~ architecture + (1 | patient), data = x, REML = FALSE), error = function(e) NULL)
    if (is.null(fit_null) || is.null(fit_full)) {
      return(tibble::tibble(n_ecDNA = n_ecDNA, n_ChrAmp = n_ChrAmp, n_patients = n_patients, chi_square = NA_real_, df = NA_real_, p_value = NA_real_, singular = NA))
    }
    lrt <- stats::anova(fit_null, fit_full)
    tibble::tibble(
      n_ecDNA = n_ecDNA,
      n_ChrAmp = n_ChrAmp,
      n_patients = n_patients,
      chi_square = lrt[2, "Chisq"],
      df = lrt[2, "Df"],
      p_value = lrt[2, "Pr(>Chisq)"],
      singular = lme4::isSingular(fit_full, tol = 1e-4)
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(FDR = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(FDR)

# -----------------------------------------------------------------------------
# 19E. Mixed-effects continuous minimum CN vs RNA
# log2(TPM + 1) ~ log2(minimum CN + 1) + (1 | patient)
# -----------------------------------------------------------------------------

mixed_cn_rna <- analysis_data |>
  dplyr::filter(!is.na(patient), patient != "", !is.na(log2_TPM), !is.na(log2_minCN)) |>
  dplyr::group_by(gene) |>
  dplyr::group_modify(\(.x, .y) {
    x <- .x
    n_patients <- dplyr::n_distinct(x$patient)
    if (nrow(x) < 4 || n_patients < 2 || dplyr::n_distinct(x$log2_minCN) < 2) {
      return(tibble::tibble(n = nrow(x), n_patients = n_patients, beta = NA_real_, chi_square = NA_real_, df = NA_real_, p_value = NA_real_, singular = NA))
    }
    fit_null <- tryCatch(lme4::lmer(log2_TPM ~ 1 + (1 | patient), data = x, REML = FALSE), error = function(e) NULL)
    fit_full <- tryCatch(lme4::lmer(log2_TPM ~ log2_minCN + (1 | patient), data = x, REML = FALSE), error = function(e) NULL)
    if (is.null(fit_null) || is.null(fit_full)) {
      return(tibble::tibble(n = nrow(x), n_patients = n_patients, beta = NA_real_, chi_square = NA_real_, df = NA_real_, p_value = NA_real_, singular = NA))
    }
    lrt <- stats::anova(fit_null, fit_full)
    tibble::tibble(
      n = nrow(x),
      n_patients = n_patients,
      beta = unname(lme4::fixef(fit_full)["log2_minCN"]),
      chi_square = lrt[2, "Chisq"],
      df = lrt[2, "Df"],
      p_value = lrt[2, "Pr(>Chisq)"],
      singular = lme4::isSingular(fit_full, tol = 1e-4)
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(FDR = stats::p.adjust(p_value, method = "BH")) |>
  dplyr::arrange(FDR)



