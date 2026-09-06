# ============================================================
# RMS WGS ecDNA analysis
# Saoud et al., CCR 2026
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(ggsci)
library(scales)
library(stringr)

broad_amplicon_colors <- c(
  "ecDNA" = "#6C77AD",
  "Chromosomal" = "#EE0000FF",
  "No focal amplification" = "#008B45FF"
)

ecdna_status_colors <- c(
  "ecDNA" = "#6C77AD",
  "non-ecDNA" = "#78c679"
)

detailed_amplicon_colors <- c(
  "ecDNA" = "#6C77AD",
  "ecDNA + FAN" = "#78c679",
  "FAN" = "#feb24c",
  "BFB" = "#CC79A7",
  "Complex-non-cyclic" = "#CE3D32",
  "Linear" = "#F0E685",
  "No focal amplification" = "#008B45FF",
  "Chromosomal" = "#EE0000FF",
  "FAN + BFB" = "#A65628",
  "ecDNA + BFB" = "#6A3D9A",
  "ecDNA + FAN + BFB" = "#8C564B"
)

result_table_path <- file.path("data", "combined_result_tables_AC2.0.0.tsv")
clinical_path <- file.path("data", "WGS_RMS_noMRN.csv")
ccdi_roster_path <- file.path("data", "samples_ccdi.xlsx")
ccdi_results_path <- file.path("data", "combined_result_tables_AC2.0.0_SJ.tsv")

data <- read.delim(result_table_path, stringsAsFactors = FALSE) %>%
  mutate(patient = sub("-T.*$", "", Sample.name))

data_clinical <- read.csv(clinical_path, stringsAsFactors = FALSE)
samples_ccdi <- readxl::read_excel(ccdi_roster_path)
ccdi_results <- read.delim(ccdi_results_path, stringsAsFactors = FALSE)

samples_ccdi <- samples_ccdi %>%
  mutate(
    CancerSubclass = case_when(
      Patient.ID %in% c("SJ000026", "SJ032029") ~ "ERMS",
      CancerSubclass == "BERMS" ~ "ERMS",
      TRUE ~ CancerSubclass
    )
  )

summary_by_disease <- data_clinical %>%
  group_by(Disease) %>%
  summarise(
    n_samples = n_distinct(DMP, na.rm = TRUE),
    n_patients = n_distinct(Individual.DMP.ID, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients))

cancer_type_summary_ccdi <- samples_ccdi %>%
  distinct(Patient.ID, Biosample.ID, CancerSubclass) %>%
  group_by(CancerSubclass) %>%
  summarise(
    n_samples = n_distinct(Biosample.ID, na.rm = TRUE),
    n_patients = n_distinct(Patient.ID, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients))

ccdi_filtered <- ccdi_results %>%
  filter(Sample.name %in% samples_ccdi$Biosample.ID) %>%
  left_join(
    samples_ccdi %>% select(Biosample.ID, Patient.ID),
    by = c("Sample.name" = "Biosample.ID")
  ) %>%
  mutate(patient = Patient.ID) %>%
  select(-Patient.ID)

joined <- bind_rows(data, ccdi_filtered)

clinical_subtype <- data_clinical %>%
  transmute(Sample.name = DMP, clinical_subtype = .data[["OncoTree.Code"]])

ccdi_subtype <- samples_ccdi %>%
  transmute(Sample.name = Biosample.ID, ccdi_subtype = CancerSubclass)

joined <- joined %>%
  left_join(ccdi_subtype, by = "Sample.name") %>%
  left_join(clinical_subtype, by = "Sample.name") %>%
  mutate(
    cancer_subtype = coalesce(ccdi_subtype, clinical_subtype),
    cancer_subtype = recode(cancer_subtype, "BERMS" = "ERMS", "SCSRMS" = "SCRMS"),
    cancer_subtype = ifelse(patient == "P-0085616" & is.na(cancer_subtype), "ARMS", cancer_subtype),
    Classification = ifelse(
      is.na(Classification) | Classification == "",
      "No focal amplification",
      Classification
    ),
    FeatureClass = case_when(
      stringr::str_detect(tidyr::replace_na(Feature.ID, ""), "^FAN_") ~ "FAN",
      stringr::str_detect(tidyr::replace_na(Feature.ID, ""), "^ecDNA_") ~ "ecDNA",
      TRUE ~ Classification
    )
  ) %>%
  select(-ccdi_subtype, -clinical_subtype)

amplicon_flags <- joined %>%
  dplyr::group_by(Sample.name, AA.amplicon.number) %>%
  dplyr::summarise(
    ecDNA_amplicon = any(FeatureClass == "ecDNA", na.rm = TRUE),
    FAN_amplicon = any(FeatureClass == "FAN", na.rm = TRUE),
    focal_amplicon = any(
      FeatureClass %in% c("ecDNA", "FAN", "BFB", "Complex-non-cyclic", "Linear"),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

joined <- joined %>%
  dplyr::left_join(
    amplicon_flags,
    by = c("Sample.name", "AA.amplicon.number")
  ) %>%
  dplyr::mutate(
    ecDNA_amplicon = tidyr::replace_na(ecDNA_amplicon, FALSE),
    FAN_amplicon = tidyr::replace_na(FAN_amplicon, FALSE),
    focal_amplicon = tidyr::replace_na(focal_amplicon, FALSE),
    ecDNA_FAN_status = dplyr::case_when(
      ecDNA_amplicon & FAN_amplicon ~ "ecDNA + FAN",
      ecDNA_amplicon ~ "ecDNA only",
      FAN_amplicon ~ "FAN only",
      TRUE ~ "Neither"
    )
  )

cohort_summary <- joined %>%
  group_by(cancer_subtype) %>%
  summarise(
    n_samples = n_distinct(Sample.name, na.rm = TRUE),
    n_patients = n_distinct(patient, na.rm = TRUE),
    .groups = "drop"
  )

classify_amplicon <- function(x) {
  case_when(
    any(x == "ecDNA", na.rm = TRUE) ~ "ecDNA",
    any(x %in% c("FAN", "BFB", "Complex-non-cyclic", "Linear"), na.rm = TRUE) ~ "Chromosomal",
    TRUE ~ "No focal amplification"
  )
}

patient_level <- joined %>%
  filter(!is.na(patient), cancer_subtype %in% c("ARMS", "ERMS")) %>%
  group_by(patient, cancer_subtype) %>%
  summarise(ecDNA_any = any(ecDNA_amplicon, na.rm = TRUE), .groups = "drop") %>%
  mutate(ecDNA_status = ifelse(ecDNA_any, "ecDNA", "non-ecDNA"))

patient_tab <- patient_level %>%
  dplyr::count(ecDNA_status, cancer_subtype) %>%
  pivot_wider(names_from = cancer_subtype, values_from = n, values_fill = 0)

fisher_result <- fisher.test(as.matrix(patient_tab[, c("ARMS", "ERMS")]))

sample_level <- joined %>%
  filter(!is.na(patient), cancer_subtype %in% c("ARMS", "ERMS")) %>%
  group_by(patient, Sample.name, cancer_subtype) %>%
  summarise(ecDNA_any = any(ecDNA_amplicon, na.rm = TRUE), .groups = "drop") %>%
  mutate(cancer_subtype = relevel(factor(cancer_subtype), ref = "ERMS"))

fit <- lme4::glmer(
  ecDNA_any ~ cancer_subtype + (1 | patient),
  family = binomial,
  data = sample_level
)

beta <- lme4::fixef(fit)["cancer_subtypeARMS"]
se <- sqrt(vcov(fit)["cancer_subtypeARMS", "cancer_subtypeARMS"])
pval <- summary(fit)$coefficients["cancer_subtypeARMS", "Pr(>|z|)"]

glmm_results <- data.frame(
  comparison = "ARMS vs ERMS",
  OR = exp(beta),
  CI_low = exp(beta - 1.96 * se),
  CI_high = exp(beta + 1.96 * se),
  p_value = pval
)

sample_level_class <- joined %>%
  filter(!is.na(Sample.name), !is.na(cancer_subtype)) %>%
  group_by(patient, Sample.name, cancer_subtype) %>%
  summarise(
    ecDNA_any = any(ecDNA_amplicon, na.rm = TRUE),
    chromosomal_any = any(
      FAN_amplicon |
        FeatureClass %in% c("FAN", "BFB", "Complex-non-cyclic", "Linear"),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Classification = case_when(
      ecDNA_any ~ "ecDNA",
      chromosomal_any ~ "Chromosomal",
      TRUE ~ "No focal amplification"
    ),
    Classification = factor(
      Classification,
      levels = c("No focal amplification", "Chromosomal", "ecDNA")
    )
  )

sample_counts <- sample_level_class %>%
  dplyr::count(cancer_subtype, Classification, name = "n")

p_sample <- ggplot(sample_counts, aes(x = cancer_subtype, y = n, fill = Classification)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(values = broad_amplicon_colors, drop = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Cancer subtype", y = "Proportion of samples", fill = "Amplicon class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

amplicon_level_class <- joined %>%
  dplyr::filter(!is.na(Sample.name), !is.na(cancer_subtype)) %>%
  dplyr::group_by(patient, Sample.name, cancer_subtype, AA.amplicon.number) %>%
  dplyr::summarise(
    ecDNA_positive = any(ecDNA_amplicon, na.rm = TRUE),
    FAN_positive = any(FAN_amplicon, na.rm = TRUE),
    chromosomal_positive = any(
      FAN_amplicon |
        FeatureClass %in% c("FAN", "BFB", "Complex-non-cyclic", "Linear"),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Classification = dplyr::case_when(
      ecDNA_positive ~ "ecDNA",
      chromosomal_positive ~ "Chromosomal",
      TRUE ~ "No focal amplification"
    ),
    ecDNA_FAN_status = dplyr::case_when(
      ecDNA_positive & FAN_positive ~ "ecDNA + FAN",
      ecDNA_positive ~ "ecDNA only",
      FAN_positive ~ "FAN only",
      TRUE ~ "Neither"
    )
  )

amplicon_counts <- amplicon_level_class %>%
  dplyr::count(cancer_subtype, Classification, name = "n") %>%
  dplyr::mutate(
    cancer_subtype = factor(cancer_subtype, levels = c("ARMS", "ERMS", "NOS", "SCRMS"))
  )

ecDNA_FAN_overlap_summary <- amplicon_level_class %>%
  dplyr::filter(ecDNA_FAN_status != "Neither") %>%
  dplyr::count(ecDNA_FAN_status, name = "n_amplicons")

cohort_amplicons <- joined %>%
  dplyr::filter(
    !is.na(Sample.name),
    !is.na(AA.amplicon.number),
    focal_amplicon
  ) %>%
  dplyr::group_by(patient, Sample.name, cancer_subtype, AA.amplicon.number) %>%
  dplyr::summarise(
    ecDNA = any(FeatureClass == "ecDNA", na.rm = TRUE),
    FAN = any(FeatureClass == "FAN", na.rm = TRUE),
    BFB = any(FeatureClass == "BFB", na.rm = TRUE),
    Complex_non_cyclic = any(FeatureClass == "Complex-non-cyclic", na.rm = TRUE),
    Linear = any(FeatureClass == "Linear", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Amplicon_classification = paste(
      c("ecDNA", "FAN", "BFB", "Complex-non-cyclic", "Linear")[
        c(ecDNA, FAN, BFB, Complex_non_cyclic, Linear)
      ],
      collapse = " + "
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(Amplicon_classification != "")

total_amplicons <- nrow(cohort_amplicons)

cohort_amplicon_classification_summary <- cohort_amplicons %>%
  dplyr::count(Amplicon_classification, name = "n_amplicons") %>%
  dplyr::mutate(
    proportion = n_amplicons / total_amplicons
  ) %>%
  dplyr::arrange(dplyr::desc(n_amplicons))

p_cohort_amplicons <- ggplot2::ggplot(
  cohort_amplicon_classification_summary,
  ggplot2::aes(
    x = reorder(Amplicon_classification, n_amplicons),
    y = n_amplicons,
    fill = Amplicon_classification
  )
) +
  ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_amplicons),
    hjust = -0.15,
    size = 4
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(values = detailed_amplicon_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    x = "Amplicon classification",
    y = "Number of amplicons",
    title = paste0("Focal amplicons in the RMS cohort (n = ", total_amplicons, ")")
  )

samples_no_focal <- joined %>%
  dplyr::filter(!is.na(Sample.name), !is.na(cancer_subtype)) %>%
  dplyr::group_by(patient, Sample.name, cancer_subtype) %>%
  dplyr::summarise(
    has_focal = any(focal_amplicon, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(!has_focal) %>%
  dplyr::transmute(
    patient,
    Sample.name,
    cancer_subtype,
    AA.amplicon.number = NA,
    Amplicon_classification = "No focal amplification"
  )

amplicons_with_no_focal <- cohort_amplicons %>%
  dplyr::select(
    patient,
    Sample.name,
    cancer_subtype,
    AA.amplicon.number,
    Amplicon_classification
  ) %>%
  dplyr::bind_rows(samples_no_focal)

amplicon_class_by_subtype <- amplicons_with_no_focal %>%
  dplyr::filter(!is.na(cancer_subtype)) %>%
  dplyr::count(
    cancer_subtype,
    Amplicon_classification,
    name = "n"
  ) %>%
  dplyr::group_by(cancer_subtype) %>%
  dplyr::mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cancer_subtype = factor(
      cancer_subtype,
      levels = c("ARMS", "ERMS", "NOS")
    ),
    Amplicon_classification = factor(
      Amplicon_classification,
      levels = c(
        "No focal amplification",
        "Linear",
        "Complex-non-cyclic",
        "BFB",
        "FAN",
        "ecDNA",
        "ecDNA + FAN",
        "FAN + BFB",
        "ecDNA + BFB",
        "ecDNA + FAN + BFB"
      )
    )
  )

p_amplicon_detailed_subtype_prop <- ggplot2::ggplot(
  amplicon_class_by_subtype,
  ggplot2::aes(
    x = cancer_subtype,
    y = n,
    fill = Amplicon_classification
  )
) +
  ggplot2::geom_col(position = "fill", width = 0.75) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(values = detailed_amplicon_colors, drop = FALSE) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    x = "Cancer subtype",
    y = "Proportion of amplicons",
    fill = "Amplicon class"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_amplicon_detailed_subtype_count <- ggplot2::ggplot(
  amplicon_class_by_subtype,
  ggplot2::aes(
    x = cancer_subtype,
    y = n,
    fill = Amplicon_classification
  )
) +
  ggplot2::geom_col(position = "stack", width = 0.75) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3
  ) +
  ggplot2::scale_fill_manual(values = detailed_amplicon_colors, drop = FALSE) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    x = "Cancer subtype",
    y = "Number of amplicons",
    fill = "Amplicon class"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

cohort_samples <- joined %>%
  dplyr::filter(!is.na(Sample.name)) %>%
  dplyr::group_by(patient, Sample.name, cancer_subtype) %>%
  dplyr::summarise(
    ecDNA = any(ecDNA_amplicon, na.rm = TRUE),
    FAN = any(FAN_amplicon, na.rm = TRUE),
    BFB = any(FeatureClass == "BFB", na.rm = TRUE),
    Complex_non_cyclic = any(FeatureClass == "Complex-non-cyclic", na.rm = TRUE),
    Linear = any(FeatureClass == "Linear", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Sample_classification = {
      cls <- c("ecDNA", "FAN", "BFB", "Complex-non-cyclic", "Linear")[
        c(ecDNA, FAN, BFB, Complex_non_cyclic, Linear)
      ]
      if (length(cls) == 0) "No focal amplification" else paste(cls, collapse = " + ")
    },
    Broad_classification = dplyr::case_when(
      ecDNA ~ "ecDNA",
      FAN | BFB | Complex_non_cyclic | Linear ~ "Chromosomal",
      TRUE ~ "No focal amplification"
    )
  ) %>%
  dplyr::ungroup()

total_samples <- dplyr::n_distinct(cohort_samples$Sample.name)

sample_level_classification_summary <- cohort_samples %>%
  dplyr::count(Sample_classification, name = "n_samples") %>%
  dplyr::mutate(proportion = n_samples / total_samples) %>%
  dplyr::arrange(dplyr::desc(n_samples))

sample_level_broad_classification_summary <- cohort_samples %>%
  dplyr::count(Broad_classification, name = "n_samples") %>%
  dplyr::mutate(
    proportion = n_samples / total_samples,
    Broad_classification = factor(
      Broad_classification,
      levels = c("No focal amplification", "Chromosomal", "ecDNA")
    )
  ) %>%
  dplyr::arrange(Broad_classification)

p_sample_detailed <- ggplot2::ggplot(
  sample_level_classification_summary,
  ggplot2::aes(
    x = reorder(Sample_classification, n_samples),
    y = n_samples,
    fill = Sample_classification
  )
) +
  ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_samples),
    hjust = -0.15,
    size = 4
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(values = detailed_amplicon_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    x = "Sample classification",
    y = "Number of samples",
    title = paste0("Sample-level amplification classes (n = ", total_samples, ")")
  )

p_sample_broad <- ggplot2::ggplot(
  sample_level_broad_classification_summary,
  ggplot2::aes(
    x = Broad_classification,
    y = n_samples,
    fill = Broad_classification
  )
) +
  ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_samples),
    vjust = -0.4,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = broad_amplicon_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    x = NULL,
    y = "Number of samples",
    title = paste0("Broad sample-level classification (n = ", total_samples, ")")
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

p_amplicons <- ggplot2::ggplot(
  amplicon_counts,
  ggplot2::aes(x = cancer_subtype, y = n, fill = Classification)
) +
  ggplot2::geom_col(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(values = broad_amplicon_colors, drop = FALSE) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    x = "Cancer subtype",
    y = "Proportion",
    fill = "Amplicon class"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

amplicon_per_sample <- joined %>%
  dplyr::filter(!is.na(Sample.name), !is.na(cancer_subtype)) %>%
  dplyr::group_by(patient, Sample.name, cancer_subtype) %>%
  dplyr::summarise(
    n_amplicons = dplyr::n_distinct(
      AA.amplicon.number[
        !is.na(AA.amplicon.number) &
          focal_amplicon
      ]
    ),
    .groups = "drop"
  )

amplicon_number_counts <- amplicon_per_sample %>%
  dplyr::count(cancer_subtype, n_amplicons, name = "n") %>%
  dplyr::mutate(
    n_amplicons = factor(
      n_amplicons,
      levels = sort(unique(as.integer(n_amplicons)))
    )
  )

p_amplicon_number <- ggplot2::ggplot(
  amplicon_number_counts,
  ggplot2::aes(x = cancer_subtype, y = n, fill = n_amplicons)
) +
  ggplot2::geom_col(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggsci::scale_fill_d3(palette = "category20c") +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::labs(
    x = "Cancer subtype",
    y = "Proportion of samples",
    fill = "Number of amplicons"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

amplicon_summary_all <- amplicon_per_sample %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean = mean(n_amplicons),
    median = median(n_amplicons),
    min = min(n_amplicons),
    max = max(n_amplicons)
  )

amplicon_summary_subtype <- amplicon_per_sample %>%
  dplyr::group_by(cancer_subtype) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean = mean(n_amplicons),
    median = median(n_amplicons),
    min = min(n_amplicons),
    max = max(n_amplicons),
    .groups = "drop"
  )

amplicon_model_data <- amplicon_per_sample %>%
  dplyr::filter(!is.na(patient), !is.na(cancer_subtype)) %>%
  dplyr::mutate(cancer_subtype = droplevels(factor(cancer_subtype)))

fit_amplicon_null <- lme4::glmer.nb(
  n_amplicons ~ 1 + (1 | patient),
  data = amplicon_model_data
)

fit_amplicon_full <- lme4::glmer.nb(
  n_amplicons ~ cancer_subtype + (1 | patient),
  data = amplicon_model_data
)

amplicon_lrt <- anova(
  fit_amplicon_null,
  fit_amplicon_full,
  test = "Chisq"
)

amplicon_burden_results <- data.frame(
  comparison = "Overall difference across RMS subtypes",
  chi_square = amplicon_lrt[2, "Chisq"],
  df = amplicon_lrt[2, "Df"],
  p_value = amplicon_lrt[2, "Pr(>Chisq)"]
)

# ============================================================
# Gene content of focal amplicons
# ============================================================

clean_list_column <- function(x) {
  x <- tidyr::replace_na(x, "")
  stringr::str_remove_all(x, "\\[|\\]|'|\"")
}

chromosomal_classes <- c("FAN", "BFB", "Complex-non-cyclic", "Linear")
focal_classes <- c("ecDNA", chromosomal_classes)

analysis_features <- joined %>%
  dplyr::mutate(
    PrimaryFeatureClass = dplyr::case_when(
      FeatureClass == "ecDNA" ~ "ecDNA",
      FeatureClass == "FAN" & !ecDNA_amplicon ~ "FAN",
      FeatureClass %in% c("BFB", "Complex-non-cyclic", "Linear") &
        !ecDNA_amplicon & !FAN_amplicon ~ FeatureClass,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(PrimaryFeatureClass))

chromosomal_gene_long <- analysis_features %>%
  dplyr::filter(PrimaryFeatureClass %in% chromosomal_classes) %>%
  dplyr::mutate(Classification = PrimaryFeatureClass) %>%
  dplyr::mutate(All.genes = clean_list_column(All.genes)) %>%
  tidyr::separate_rows(All.genes, sep = ",") %>%
  dplyr::mutate(gene = stringr::str_trim(All.genes)) %>%
  dplyr::filter(gene != "") %>%
  dplyr::distinct(patient, Sample.name, Classification, gene)

recurrent_chromosomal_genes <- chromosomal_gene_long %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_samples = dplyr::n_distinct(Sample.name),
    classifications = paste(sort(unique(Classification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_samples > 1) %>%
  dplyr::arrange(dplyr::desc(n_samples), gene)

recurrent_chromosomal_genes_patient <- chromosomal_gene_long %>%
  dplyr::filter(!is.na(patient)) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_patients = dplyr::n_distinct(patient),
    classifications = paste(sort(unique(Classification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_patients > 1) %>%
  dplyr::arrange(dplyr::desc(n_patients), gene)

chromosomal_oncogene_long <- analysis_features %>%
  dplyr::filter(PrimaryFeatureClass %in% chromosomal_classes) %>%
  dplyr::mutate(Classification = PrimaryFeatureClass) %>%
  dplyr::mutate(Oncogenes = clean_list_column(Oncogenes)) %>%
  tidyr::separate_rows(Oncogenes, sep = ",") %>%
  dplyr::mutate(oncogene = stringr::str_trim(Oncogenes)) %>%
  dplyr::filter(oncogene != "") %>%
  dplyr::distinct(patient, Sample.name, Classification, oncogene)

recurrent_chromosomal_oncogenes <- chromosomal_oncogene_long %>%
  dplyr::group_by(oncogene) %>%
  dplyr::summarise(
    n_samples = dplyr::n_distinct(Sample.name),
    classifications = paste(sort(unique(Classification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_samples > 1) %>%
  dplyr::arrange(dplyr::desc(n_samples), oncogene)

recurrent_chromosomal_oncogenes_patient <- chromosomal_oncogene_long %>%
  dplyr::filter(!is.na(patient)) %>%
  dplyr::group_by(oncogene) %>%
  dplyr::summarise(
    n_patients = dplyr::n_distinct(patient),
    classifications = paste(sort(unique(Classification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_patients > 1) %>%
  dplyr::arrange(dplyr::desc(n_patients), oncogene)

# ============================================================
# Exploratory gene enrichment in ecDNA versus chromosomal amplicon features
# ============================================================

focal_features <- analysis_features %>%
  dplyr::filter(
    PrimaryFeatureClass %in% focal_classes,
    !is.na(Sample.name),
    !is.na(Feature.ID)
  ) %>%
  dplyr::mutate(
    Classification = PrimaryFeatureClass,
    AmpliconClass = ifelse(PrimaryFeatureClass == "ecDNA", "ecDNA", "Chromosomal")
  ) %>%
  dplyr::distinct(patient, Sample.name, Feature.ID, AmpliconClass, .keep_all = TRUE)

feature_gene_long <- focal_features %>%
  dplyr::mutate(All.genes = clean_list_column(All.genes)) %>%
  tidyr::separate_rows(All.genes, sep = ",") %>%
  dplyr::mutate(gene = stringr::str_trim(All.genes)) %>%
  dplyr::filter(gene != "") %>%
  dplyr::distinct(Sample.name, Feature.ID, AmpliconClass, gene)

feature_totals <- focal_features %>%
  dplyr::count(AmpliconClass, name = "n_features")

n_ecDNA_features <- feature_totals$n_features[feature_totals$AmpliconClass == "ecDNA"]
n_chromosomal_features <- feature_totals$n_features[feature_totals$AmpliconClass == "Chromosomal"]

gene_assoc <- feature_gene_long %>%
  dplyr::count(gene, AmpliconClass, name = "present") %>%
  tidyr::complete(
    gene,
    AmpliconClass = c("ecDNA", "Chromosomal"),
    fill = list(present = 0)
  ) %>%
  tidyr::pivot_wider(names_from = AmpliconClass, values_from = present) %>%
  dplyr::filter(ecDNA + Chromosomal >= 2) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    test = list(stats::fisher.test(matrix(
      c(
        ecDNA, n_ecDNA_features - ecDNA,
        Chromosomal, n_chromosomal_features - Chromosomal
      ),
      nrow = 2,
      byrow = TRUE
    ))),
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR = stats::p.adjust(p_value, method = "BH"),
    log2OR = log2(odds_ratio)
  ) %>%
  dplyr::arrange(FDR)

top_genes <- gene_assoc %>%
  dplyr::slice_head(n = 40) %>%
  dplyr::pull(gene)

gene_plot_data <- gene_assoc %>%
  dplyr::filter(gene %in% top_genes) %>%
  tidyr::pivot_longer(
    cols = c(ecDNA, Chromosomal),
    names_to = "AmpliconClass",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    log2OR_plot = pmax(pmin(log2OR, 3), -3),
    color_value = ifelse(AmpliconClass == "ecDNA", log2OR_plot, NA_real_)
  )

p_gene <- ggplot2::ggplot(
  gene_plot_data,
  ggplot2::aes(x = AmpliconClass, y = reorder(gene, log2OR))
) +
  ggplot2::geom_point(
    ggplot2::aes(size = ifelse(count == 0, NA, count), color = color_value)
  ) +
  ggplot2::scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    na.value = "grey80"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    color = "ecDNA enrichment\n(log2 OR)",
    size = "Feature count"
  )

# ============================================================
# ecDNA features without an annotated oncogene
# ============================================================

has_oncogene <- function(x) {
  x <- stringr::str_trim(clean_list_column(x))
  !(x %in% c("", "NA", "NULL"))
}

ecdna_no_oncogene_genes <- analysis_features %>%
  dplyr::filter(PrimaryFeatureClass == "ecDNA") %>%
  dplyr::mutate(
    oncogene_present = has_oncogene(Oncogenes),
    All.genes = clean_list_column(All.genes)
  ) %>%
  dplyr::filter(!oncogene_present) %>%
  tidyr::separate_rows(All.genes, sep = ",") %>%
  dplyr::mutate(gene = stringr::str_trim(All.genes)) %>%
  dplyr::filter(gene != "") %>%
  dplyr::distinct(Sample.name, gene) %>%
  dplyr::count(gene, name = "n_samples", sort = TRUE)

recurrent_ecdna_no_oncogene_genes <- ecdna_no_oncogene_genes %>%
  dplyr::filter(n_samples > 1)

# ============================================================
# Amplicon complexity
# ============================================================

complexity_df <- analysis_features %>%
  dplyr::filter(
    PrimaryFeatureClass %in% focal_classes,
    !is.na(patient),
    !is.na(Feature.ID)
  ) %>%
  dplyr::mutate(
    Complexity.score = as.numeric(Complexity.score),
    Classification = factor(
      PrimaryFeatureClass,
      levels = c("BFB", "Complex-non-cyclic", "FAN", "ecDNA", "Linear")
    )
  ) %>%
  dplyr::filter(is.finite(Complexity.score), !is.na(Classification)) %>%
  dplyr::distinct(patient, Sample.name, Feature.ID, Classification, .keep_all = TRUE)

complexity_summary <- complexity_df %>%
  dplyr::group_by(Classification) %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    median_complexity = median(Complexity.score),
    .groups = "drop"
  )

fit_complexity_null <- lme4::lmer(
  log1p(Complexity.score) ~ 1 + (1 | patient),
  data = complexity_df,
  REML = FALSE
)

fit_complexity_full <- lme4::lmer(
  log1p(Complexity.score) ~ Classification + (1 | patient),
  data = complexity_df,
  REML = FALSE
)

complexity_lrt <- anova(fit_complexity_null, fit_complexity_full)

complexity_test <- data.frame(
  comparison = "Overall difference across amplicon classes",
  chi_square = complexity_lrt[2, "Chisq"],
  df = complexity_lrt[2, "Df"],
  p_value = complexity_lrt[2, "Pr(>Chisq)"]
)

complexity_p_label <- ifelse(
  complexity_test$p_value < 0.001,
  "Global P < 0.001",
  paste0("Global P = ", formatC(complexity_test$p_value, format = "f", digits = 3))
)

complexity_levels <- levels(droplevels(complexity_df$Classification))
complexity_pairs <- combn(complexity_levels, 2, simplify = FALSE)

complexity_pairwise <- lapply(complexity_pairs, function(pair) {
  pair_df <- complexity_df %>%
    dplyr::filter(Classification %in% pair) %>%
    dplyr::mutate(Classification = droplevels(Classification))

  fit_null <- lme4::lmer(
    log1p(Complexity.score) ~ 1 + (1 | patient),
    data = pair_df,
    REML = FALSE
  )

  fit_full <- lme4::lmer(
    log1p(Complexity.score) ~ Classification + (1 | patient),
    data = pair_df,
    REML = FALSE
  )

  lrt <- anova(fit_null, fit_full)

  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    p_value = lrt[2, "Pr(>Chisq)"]
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    significance = dplyr::case_when(
      FDR < 0.0001 ~ "****",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x1 = match(group1, complexity_levels),
    x2 = match(group2, complexity_levels)
  ) %>%
  dplyr::arrange(x2 - x1, x1)

complexity_range <- diff(range(complexity_df$Complexity.score, na.rm = TRUE))
if (complexity_range == 0) complexity_range <- 1

complexity_pairwise <- complexity_pairwise %>%
  dplyr::mutate(
    y = max(complexity_df$Complexity.score, na.rm = TRUE) +
      (0.16 + (dplyr::row_number() - 1) * 0.11) * complexity_range,
    tip = 0.025 * complexity_range
  )

global_y <- max(complexity_df$Complexity.score, na.rm = TRUE) + 0.05 * complexity_range
upper_y <- max(complexity_pairwise$y, na.rm = TRUE) + 0.10 * complexity_range

p_complexity <- ggplot2::ggplot(
  complexity_df,
  ggplot2::aes(x = Classification, y = Complexity.score, fill = Classification)
) +
  ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.80) +
  ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  ggplot2::geom_jitter(width = 0.10, alpha = 0.25, size = 0.9) +
  ggplot2::geom_segment(
    data = complexity_pairwise,
    ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = complexity_pairwise,
    ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tip),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = complexity_pairwise,
    ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tip),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = complexity_pairwise,
    ggplot2::aes(x = (x1 + x2) / 2, y = y + 0.02 * complexity_range, label = significance),
    inherit.aes = FALSE,
    size = 3.8,
    vjust = 0
  ) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = global_y,
    label = complexity_p_label,
    hjust = 0,
    vjust = 0,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = detailed_amplicon_colors, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(NA, upper_y), clip = "off") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(x = NULL, y = "Complexity score", fill = "Amplicon class") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.margin = ggplot2::margin(10, 15, 10, 10)
  )

# ============================================================
# Oncogene content by parent amplicon class
# ============================================================

parent_oncogene_status <- joined %>%
  dplyr::filter(
    !is.na(patient),
    !is.na(Sample.name),
    !is.na(AA.amplicon.number),
    focal_amplicon
  ) %>%
  dplyr::group_by(patient, Sample.name, AA.amplicon.number) %>%
  dplyr::summarise(
    oncogene_present = any(has_oncogene(Oncogenes), na.rm = TRUE),
    .groups = "drop"
  )

oncogene_df <- cohort_amplicons %>%
  dplyr::left_join(
    parent_oncogene_status,
    by = c("patient", "Sample.name", "AA.amplicon.number")
  ) %>%
  dplyr::mutate(
    oncogene_present = tidyr::replace_na(oncogene_present, FALSE),
    Classification = dplyr::case_when(
      ecDNA ~ "ecDNA",
      FAN ~ "FAN",
      BFB ~ "BFB",
      Complex_non_cyclic ~ "Complex-non-cyclic",
      Linear ~ "Linear",
      TRUE ~ NA_character_
    ),
    Classification = factor(
      Classification,
      levels = c("ecDNA", "FAN", "BFB", "Complex-non-cyclic", "Linear")
    )
  ) %>%
  dplyr::filter(!is.na(Classification))

oncogene_counts <- oncogene_df %>%
  dplyr::count(Classification, oncogene_present, name = "n") %>%
  dplyr::mutate(
    oncogene_group = factor(
      oncogene_present,
      levels = c(FALSE, TRUE),
      labels = c("No oncogene", "Oncogene")
    )
  )

p_oncogene <- ggplot2::ggplot(
  oncogene_counts,
  ggplot2::aes(x = Classification, y = n, fill = oncogene_group)
) +
  ggplot2::geom_col(position = "fill", width = 0.7) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::scale_fill_manual(values = c("No oncogene" = "grey80", "Oncogene" = "#D55E00")) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(x = NULL, y = "Proportion of amplicons", fill = NULL) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

p_oncogene_count <- ggplot2::ggplot(
  oncogene_counts,
  ggplot2::aes(x = Classification, y = n, fill = oncogene_group)
) +
  ggplot2::geom_col(position = "stack", width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_stack(vjust = 0.5),
    size = 3.2
  ) +
  ggplot2::scale_fill_manual(values = c("No oncogene" = "grey80", "Oncogene" = "#D55E00")) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(x = NULL, y = "Number of amplicons", fill = NULL) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

oncogene_model_data <- oncogene_df %>%
  dplyr::filter(!is.na(patient), !is.na(oncogene_present)) %>%
  dplyr::mutate(
    Classification = droplevels(Classification),
    oncogene_present = as.integer(oncogene_present)
  )

fit_oncogene_null <- lme4::glmer(
  oncogene_present ~ 1 + (1 | patient),
  family = binomial,
  data = oncogene_model_data
)

fit_oncogene_full <- lme4::glmer(
  oncogene_present ~ Classification + (1 | patient),
  family = binomial,
  data = oncogene_model_data
)

oncogene_global_lrt <- anova(
  fit_oncogene_null,
  fit_oncogene_full,
  test = "Chisq"
)

oncogene_global_test <- data.frame(
  comparison = "Overall oncogene enrichment across amplicon classes",
  chi_square = oncogene_global_lrt[2, "Chisq"],
  df = oncogene_global_lrt[2, "Df"],
  p_value = oncogene_global_lrt[2, "Pr(>Chisq)"]
)

oncogene_levels <- levels(oncogene_model_data$Classification)
oncogene_pairs <- combn(oncogene_levels, 2, simplify = FALSE)

oncogene_pairwise_tests <- lapply(oncogene_pairs, function(pair) {
  pair_df <- oncogene_model_data %>%
    dplyr::filter(Classification %in% pair) %>%
    dplyr::mutate(Classification = droplevels(Classification))

  fit_null <- lme4::glmer(
    oncogene_present ~ 1 + (1 | patient),
    family = binomial,
    data = pair_df
  )

  fit_full <- lme4::glmer(
    oncogene_present ~ Classification + (1 | patient),
    family = binomial,
    data = pair_df
  )

  lrt <- anova(fit_null, fit_full, test = "Chisq")
  beta <- lme4::fixef(fit_full)[2]
  se <- sqrt(vcov(fit_full)[2, 2])

  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    OR_group2_vs_group1 = exp(beta),
    CI_low = exp(beta - 1.96 * se),
    CI_high = exp(beta + 1.96 * se),
    p_value = lrt[2, "Pr(>Chisq)"]
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    significance = dplyr::case_when(
      FDR < 0.0001 ~ "****",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# ============================================================
# Complexity score by oncogene presence: parent-amplicon level
# ============================================================

complexity_oncogene_df <- joined %>%
  dplyr::select(
    patient,
    Sample.name,
    AA.amplicon.number,
    FeatureClass,
    Complexity.score
  ) %>%
  dplyr::inner_join(
    oncogene_df %>%
      dplyr::select(
        patient,
        Sample.name,
        AA.amplicon.number,
        Classification,
        oncogene_present
      ),
    by = c("patient", "Sample.name", "AA.amplicon.number")
  ) %>%
  dplyr::mutate(
    Complexity.score = as.numeric(Complexity.score),
    Classification_chr = as.character(Classification)
  ) %>%
  dplyr::filter(
    FeatureClass == Classification_chr,
    is.finite(Complexity.score)
  ) %>%
  dplyr::group_by(
    patient,
    Sample.name,
    AA.amplicon.number,
    Classification,
    oncogene_present
  ) %>%
  dplyr::summarise(
    Complexity.score = median(Complexity.score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    oncogene_group = factor(
      oncogene_present,
      levels = c(FALSE, TRUE),
      labels = c("No oncogene", "Oncogene present")
    )
  )

fit_oncogene_complexity_null <- lme4::lmer(
  log1p(Complexity.score) ~ 1 + (1 | patient),
  data = complexity_oncogene_df,
  REML = FALSE
)

fit_oncogene_complexity_full <- lme4::lmer(
  log1p(Complexity.score) ~ oncogene_present + (1 | patient),
  data = complexity_oncogene_df,
  REML = FALSE
)

oncogene_complexity_lrt <- anova(
  fit_oncogene_complexity_null,
  fit_oncogene_complexity_full
)

oncogene_complexity_test <- data.frame(
  comparison = "Oncogene present vs absent",
  chi_square = oncogene_complexity_lrt[2, "Chisq"],
  df = oncogene_complexity_lrt[2, "Df"],
  p_value = oncogene_complexity_lrt[2, "Pr(>Chisq)"]
)

oncogene_complexity_p_label <- ifelse(
  oncogene_complexity_test$p_value < 0.001,
  "LRT P < 0.001",
  paste0(
    "LRT P = ",
    formatC(oncogene_complexity_test$p_value, format = "f", digits = 3)
  )
)

oncogene_complexity_signif <- dplyr::case_when(
  oncogene_complexity_test$p_value < 0.0001 ~ "****",
  oncogene_complexity_test$p_value < 0.001 ~ "***",
  oncogene_complexity_test$p_value < 0.01 ~ "**",
  oncogene_complexity_test$p_value < 0.05 ~ "*",
  TRUE ~ "ns"
)

oncogene_range <- diff(range(complexity_oncogene_df$Complexity.score, na.rm = TRUE))
if (oncogene_range == 0) oncogene_range <- 1

oncogene_global_y <- max(complexity_oncogene_df$Complexity.score, na.rm = TRUE) +
  0.05 * oncogene_range
oncogene_bracket_y <- max(complexity_oncogene_df$Complexity.score, na.rm = TRUE) +
  0.20 * oncogene_range
oncogene_tip <- 0.025 * oncogene_range
oncogene_upper_y <- oncogene_bracket_y + 0.10 * oncogene_range

p_oncogene_complexity <- ggplot2::ggplot(
  complexity_oncogene_df,
  ggplot2::aes(x = oncogene_group, y = Complexity.score, fill = oncogene_group)
) +
  ggplot2::geom_violin(trim = FALSE, scale = "width", alpha = 0.80) +
  ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  ggplot2::geom_jitter(width = 0.10, alpha = 0.25, size = 0.9) +
  ggplot2::annotate(
    "segment",
    x = 1, xend = 2,
    y = oncogene_bracket_y, yend = oncogene_bracket_y,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "segment",
    x = 1, xend = 1,
    y = oncogene_bracket_y, yend = oncogene_bracket_y - oncogene_tip,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "segment",
    x = 2, xend = 2,
    y = oncogene_bracket_y, yend = oncogene_bracket_y - oncogene_tip,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = oncogene_bracket_y + 0.02 * oncogene_range,
    label = oncogene_complexity_signif,
    size = 3.8,
    vjust = 0
  ) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = oncogene_global_y,
    label = oncogene_complexity_p_label,
    hjust = 0,
    vjust = 0,
    size = 4
  ) +
  ggsci::scale_fill_bmj() +
  ggplot2::coord_cartesian(ylim = c(NA, oncogene_upper_y), clip = "off") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(x = NULL, y = "Complexity score", fill = NULL) +
  ggplot2::theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(10, 15, 10, 10)
  )

# ecDNA context of recurrent genes
# ============================================================

genes_of_interest <- c("MYC", "MYCN", "MDM2", "CDK4", "FOXO1", "PAX7", "FGFR1")

ecdna_context <- analysis_features %>%
  dplyr::filter(
    PrimaryFeatureClass == "ecDNA",
    !is.na(ecDNA.context),
    ecDNA.context != ""
  ) %>%
  dplyr::mutate(All.genes = clean_list_column(All.genes)) %>%
  tidyr::separate_rows(All.genes, sep = ",") %>%
  dplyr::mutate(gene = stringr::str_trim(All.genes)) %>%
  dplyr::filter(gene %in% genes_of_interest) %>%
  dplyr::distinct(patient, Sample.name, Feature.ID, gene, ecDNA.context)

ecdna_context_counts <- ecdna_context %>%
  dplyr::count(gene, ecDNA.context, name = "n")

p_context <- ggplot2::ggplot(
  ecdna_context_counts,
  ggplot2::aes(x = gene, y = n, fill = ecDNA.context)
) +
  ggplot2::geom_col(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggsci::scale_fill_npg() +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(x = "Gene", y = "Proportion", fill = "ecDNA context") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

# ============================================================
# Copy number: ecDNA vs chromosomal amplification features
# ============================================================

analysis_features_no_ecDNA_FAN <- analysis_features %>%
  dplyr::filter(!(ecDNA_amplicon & FAN_amplicon))

copy_number_features <- analysis_features %>%
  dplyr::filter(
    PrimaryFeatureClass %in% focal_classes,
    !is.na(patient),
    !is.na(Feature.ID)
  ) %>%
  dplyr::mutate(
    Classification = PrimaryFeatureClass,
    Class2 = ifelse(PrimaryFeatureClass == "ecDNA", "ecDNA", "Chromosomal"),
    Feature.median.copy.number = as.numeric(Feature.median.copy.number),
    Feature.maximum.copy.number = as.numeric(Feature.maximum.copy.number)
  ) %>%
  dplyr::distinct(
    patient,
    Sample.name,
    AA.amplicon.number,
    Feature.ID,
    Class2,
    .keep_all = TRUE
  )

copy_number_summary <- copy_number_features %>%
  dplyr::group_by(Class2) %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    median_CN = median(Feature.median.copy.number, na.rm = TRUE),
    median_max_CN = median(Feature.maximum.copy.number, na.rm = TRUE),
    .groups = "drop"
  )

copy_number_model_data <- copy_number_features %>%
  dplyr::filter(
    is.finite(Feature.median.copy.number),
    Feature.median.copy.number > 0
  )

fit_cn_null <- lme4::lmer(
  log10(Feature.median.copy.number) ~ 1 + (1 | patient),
  data = copy_number_model_data,
  REML = FALSE
)

fit_cn_full <- lme4::lmer(
  log10(Feature.median.copy.number) ~ Class2 + (1 | patient),
  data = copy_number_model_data,
  REML = FALSE
)

cn_lrt <- anova(fit_cn_null, fit_cn_full)

copy_number_test <- data.frame(
  comparison = "ecDNA vs Chromosomal",
  chi_square = cn_lrt[2, "Chisq"],
  df = cn_lrt[2, "Df"],
  p_value = cn_lrt[2, "Pr(>Chisq)"]
)

copy_number_p_label <- ifelse(
  copy_number_test$p_value < 0.001,
  "P < 0.001",
  paste0("P = ", formatC(copy_number_test$p_value, format = "f", digits = 3))
)

copy_number_max <- max(copy_number_model_data$Feature.median.copy.number, na.rm = TRUE)
copy_number_bracket_y <- copy_number_max * 1.25
copy_number_tip_y <- copy_number_bracket_y / 1.06
copy_number_label_y <- copy_number_bracket_y * 1.07

class_colors <- broad_amplicon_colors[c("ecDNA", "Chromosomal")]

p_cn <- ggplot2::ggplot(
  copy_number_model_data,
  ggplot2::aes(x = Class2, y = Feature.median.copy.number, fill = Class2)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = Class2),
    width = 0.15,
    size = 1.8,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "segment",
    x = 1, xend = 2,
    y = copy_number_bracket_y, yend = copy_number_bracket_y,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "segment",
    x = 1, xend = 1,
    y = copy_number_bracket_y, yend = copy_number_tip_y,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "segment",
    x = 2, xend = 2,
    y = copy_number_bracket_y, yend = copy_number_tip_y,
    linewidth = 0.35
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = copy_number_label_y,
    label = copy_number_p_label,
    size = 4
  ) +
  ggplot2::scale_y_log10(
    expand = ggplot2::expansion(mult = c(0.05, 0.18))
  ) +
  ggplot2::scale_fill_manual(values = class_colors) +
  ggplot2::scale_color_manual(values = class_colors) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = NULL,
    y = "Median copy number (log10)"
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(10, 15, 10, 10)
  )

# ============================================================
# Copy number across detailed amplicon classes
# ============================================================

copy_number_class_data <- analysis_features %>%
  dplyr::filter(
    PrimaryFeatureClass %in% c("BFB", "Complex-non-cyclic", "FAN", "ecDNA", "Linear"),
    !is.na(patient),
    !is.na(Feature.ID)
  ) %>%
  dplyr::mutate(
    Feature.median.copy.number = as.numeric(Feature.median.copy.number),
    AmpliconClass = factor(
      PrimaryFeatureClass,
      levels = c("BFB", "Complex-non-cyclic", "FAN", "ecDNA", "Linear")
    )
  ) %>%
  dplyr::filter(
    is.finite(Feature.median.copy.number),
    Feature.median.copy.number > 0,
    !is.na(AmpliconClass)
  ) %>%
  dplyr::distinct(
    patient,
    Sample.name,
    AA.amplicon.number,
    Feature.ID,
    AmpliconClass,
    .keep_all = TRUE
  )

copy_number_class_summary <- copy_number_class_data %>%
  dplyr::group_by(AmpliconClass) %>%
  dplyr::summarise(
    n_features = dplyr::n(),
    median_CN = median(Feature.median.copy.number),
    .groups = "drop"
  )

fit_cn_class_null <- lme4::lmer(
  log10(Feature.median.copy.number) ~ 1 + (1 | patient),
  data = copy_number_class_data,
  REML = FALSE
)

fit_cn_class_full <- lme4::lmer(
  log10(Feature.median.copy.number) ~ AmpliconClass + (1 | patient),
  data = copy_number_class_data,
  REML = FALSE
)

cn_class_lrt <- anova(fit_cn_class_null, fit_cn_class_full)

copy_number_class_global_test <- data.frame(
  comparison = "Overall difference across amplicon classes",
  chi_square = cn_class_lrt[2, "Chisq"],
  df = cn_class_lrt[2, "Df"],
  p_value = cn_class_lrt[2, "Pr(>Chisq)"]
)

cn_global_label <- ifelse(
  copy_number_class_global_test$p_value < 0.001,
  "Global P < 0.001",
  paste0(
    "Global P = ",
    formatC(copy_number_class_global_test$p_value, format = "f", digits = 3)
  )
)

cn_levels <- levels(droplevels(copy_number_class_data$AmpliconClass))
cn_pairs <- combn(cn_levels, 2, simplify = FALSE)

copy_number_class_pairwise <- lapply(cn_pairs, function(pair) {
  pair_df <- copy_number_class_data %>%
    dplyr::filter(AmpliconClass %in% pair) %>%
    dplyr::mutate(AmpliconClass = droplevels(AmpliconClass))

  fit_null <- lme4::lmer(
    log10(Feature.median.copy.number) ~ 1 + (1 | patient),
    data = pair_df,
    REML = FALSE
  )

  fit_full <- lme4::lmer(
    log10(Feature.median.copy.number) ~ AmpliconClass + (1 | patient),
    data = pair_df,
    REML = FALSE
  )

  lrt <- anova(fit_null, fit_full)

  data.frame(
    group1 = pair[1],
    group2 = pair[2],
    p_value = lrt[2, "Pr(>Chisq)"]
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    significance = dplyr::case_when(
      FDR < 0.0001 ~ "****",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x1 = match(group1, cn_levels),
    x2 = match(group2, cn_levels)
  ) %>%
  dplyr::arrange(x2 - x1, x1)

cn_log_range <- diff(range(log10(copy_number_class_data$Feature.median.copy.number), na.rm = TRUE))
if (cn_log_range == 0) cn_log_range <- 1

cn_max_log <- max(log10(copy_number_class_data$Feature.median.copy.number), na.rm = TRUE)

copy_number_class_pairwise <- copy_number_class_pairwise %>%
  dplyr::mutate(
    y_log = cn_max_log + 0.12 + (dplyr::row_number() - 1) * 0.10 * cn_log_range,
    y = 10^y_log,
    tip = 10^(y_log - 0.025 * cn_log_range)
  )

cn_global_y <- 10^(cn_max_log + 0.04 * cn_log_range)
cn_upper_y <- 10^(max(copy_number_class_pairwise$y_log) + 0.08 * cn_log_range)

p_cn_class <- ggplot2::ggplot(
  copy_number_class_data,
  ggplot2::aes(
    x = AmpliconClass,
    y = Feature.median.copy.number,
    fill = AmpliconClass
  )
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  ggplot2::geom_jitter(
    width = 0.10,
    alpha = 0.35,
    size = 1
  ) +
  ggplot2::geom_segment(
    data = copy_number_class_pairwise,
    ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = copy_number_class_pairwise,
    ggplot2::aes(x = x1, xend = x1, y = y, yend = tip),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = copy_number_class_pairwise,
    ggplot2::aes(x = x2, xend = x2, y = y, yend = tip),
    inherit.aes = FALSE,
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
    data = copy_number_class_pairwise,
    ggplot2::aes(
      x = (x1 + x2) / 2,
      y = y,
      label = significance
    ),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3.6
  ) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = cn_global_y,
    label = cn_global_label,
    hjust = 0,
    vjust = 0,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = detailed_amplicon_colors, drop = FALSE) +
  ggplot2::scale_y_log10() +
  ggplot2::coord_cartesian(
    ylim = c(NA, cn_upper_y),
    clip = "off"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = NULL,
    y = "Median copy number (log10)",
    fill = "Amplicon class"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.margin = ggplot2::margin(10, 15, 10, 10)
  )

# ============================================================
# Amplification class by tumor stage and therapy timepoint
# ============================================================

required_msk_meta <- c("DMP", "Tumor.Class", "Therapy.Timepoint")
required_ccdi_meta <- c("Biosample.ID", "Tumor_history")

stopifnot(all(required_msk_meta %in% names(data_clinical)))
stopifnot(all(required_ccdi_meta %in% names(samples_ccdi)))

msk_stage_meta <- data_clinical %>%
  dplyr::distinct(DMP, .keep_all = TRUE) %>%
  dplyr::transmute(
    sample_id = DMP,
    tumor_status = as.character(Tumor.Class),
    therapy_timepoint = as.character(Therapy.Timepoint),
    source = "MSK"
  )

ccdi_stage_meta <- samples_ccdi %>%
  dplyr::distinct(Biosample.ID, .keep_all = TRUE) %>%
  dplyr::transmute(
    sample_id = Biosample.ID,
    tumor_status = as.character(Tumor_history),
    therapy_timepoint = NA_character_,
    source = "CCDI"
  )

sample_stage_data <- cohort_samples %>%
  dplyr::transmute(
    patient,
    sample_id = Sample.name,
    cancer_subtype,
    amplicon_class = as.character(Broad_classification)
  ) %>%
  dplyr::left_join(
    dplyr::bind_rows(msk_stage_meta, ccdi_stage_meta),
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    tumor_status_std = toupper(trimws(tumor_status)),
    tumor_stage = dplyr::case_when(
      tumor_status_std %in% c("PRIMARY", "DIAGNOSIS") ~ "Primary",
      !is.na(tumor_status_std) & tumor_status_std != "" ~ "Advanced",
      TRUE ~ NA_character_
    ),
    tumor_stage = factor(tumor_stage, levels = c("Primary", "Advanced")),
    amplicon_class = factor(
      amplicon_class,
      levels = c("ecDNA", "Chromosomal", "No focal amplification")
    ),
    ecDNA_present = as.integer(amplicon_class == "ecDNA")
  )

# ------------------------------------------------------------
# Primary vs advanced: all three amplification classes
# Fisher exact test is shown for the overall 2 x 3 distribution.
# ------------------------------------------------------------

stage_class_data <- sample_stage_data %>%
  dplyr::filter(!is.na(tumor_stage), !is.na(amplicon_class))

tab_stage_class <- table(
  stage_class_data$tumor_stage,
  stage_class_data$amplicon_class
)

fisher_stage_class <- fisher.test(tab_stage_class)

stage_class_test <- data.frame(
  comparison = "Primary vs Advanced: overall amplicon-class distribution",
  test = "Fisher exact",
  p_value = fisher_stage_class$p.value
)

stage_class_label <- ifelse(
  fisher_stage_class$p.value < 0.001,
  "Fisher P < 0.001",
  paste0("Fisher P = ", signif(fisher_stage_class$p.value, 3))
)

p_stage_class <- ggplot2::ggplot(
  stage_class_data,
  ggplot2::aes(x = tumor_stage, fill = amplicon_class)
) +
  ggplot2::geom_bar(position = "fill", width = 0.65) +
  ggplot2::scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.10),
    expand = c(0, 0)
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = 1.05,
    label = stage_class_label,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = broad_amplicon_colors, drop = FALSE) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion of samples",
    fill = "Amplicon class",
    title = "Amplicon class distribution in primary vs advanced tumors"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

# ------------------------------------------------------------
# Primary vs advanced: ecDNA vs non-ecDNA
# OR > 1 means higher odds of ecDNA in advanced tumors.
# ------------------------------------------------------------

stage_ecdna_data <- stage_class_data %>%
  dplyr::filter(!is.na(patient)) %>%
  dplyr::mutate(
    tumor_stage = stats::relevel(tumor_stage, ref = "Primary")
  )

fit_stage_ecdna <- lme4::glmer(
  ecDNA_present ~ tumor_stage + (1 | patient),
  family = binomial,
  data = stage_ecdna_data
)

beta_stage <- lme4::fixef(fit_stage_ecdna)["tumor_stageAdvanced"]
se_stage <- sqrt(vcov(fit_stage_ecdna)["tumor_stageAdvanced", "tumor_stageAdvanced"])
p_stage <- summary(fit_stage_ecdna)$coefficients[
  "tumor_stageAdvanced",
  "Pr(>|z|)"
]

stage_ecdna_test <- data.frame(
  comparison = "Advanced vs Primary",
  OR = exp(beta_stage),
  CI_low = exp(beta_stage - 1.96 * se_stage),
  CI_high = exp(beta_stage + 1.96 * se_stage),
  p_value = p_stage
)

stage_ecdna_label <- paste0(
  "OR = ", round(stage_ecdna_test$OR, 2),
  "\nP = ", signif(stage_ecdna_test$p_value, 3)
)

stage_ecdna_plot_data <- stage_ecdna_data %>%
  dplyr::mutate(
    ecDNA_status = factor(
      ifelse(ecDNA_present == 1, "ecDNA", "non-ecDNA"),
      levels = c("non-ecDNA", "ecDNA")
    )
  )

p_stage_ecdna <- ggplot2::ggplot(
  stage_ecdna_plot_data,
  ggplot2::aes(x = tumor_stage, fill = ecDNA_status)
) +
  ggplot2::geom_bar(position = "fill", width = 0.60) +
  ggplot2::scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.10),
    expand = c(0, 0)
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = 1.05,
    label = stage_ecdna_label,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = ecdna_status_colors, drop = FALSE) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion of samples",
    fill = "Amplicon status",
    title = "ecDNA frequency in primary vs advanced tumors"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

# ------------------------------------------------------------
# Pre-therapy vs post-therapy: MSK samples with known timepoint
# ------------------------------------------------------------

therapy_data <- sample_stage_data %>%
  dplyr::filter(
    source == "MSK",
    therapy_timepoint %in% c("PRE-THERAPY", "POST-THERAPY"),
    !is.na(amplicon_class)
  ) %>%
  dplyr::mutate(
    therapy_timepoint = factor(
      therapy_timepoint,
      levels = c("PRE-THERAPY", "POST-THERAPY")
    )
  )

tab_therapy_class <- table(
  therapy_data$therapy_timepoint,
  therapy_data$amplicon_class
)

fisher_therapy_class <- fisher.test(tab_therapy_class)

therapy_class_test <- data.frame(
  comparison = "Pre-therapy vs Post-therapy: overall amplicon-class distribution",
  test = "Fisher exact",
  p_value = fisher_therapy_class$p.value
)

therapy_class_label <- ifelse(
  fisher_therapy_class$p.value < 0.001,
  "Fisher P < 0.001",
  paste0("Fisher P = ", signif(fisher_therapy_class$p.value, 3))
)

p_therapy_class <- ggplot2::ggplot(
  therapy_data,
  ggplot2::aes(x = therapy_timepoint, fill = amplicon_class)
) +
  ggplot2::geom_bar(position = "fill", width = 0.65) +
  ggplot2::scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.10),
    expand = c(0, 0)
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = 1.05,
    label = therapy_class_label,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = broad_amplicon_colors, drop = FALSE) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion of samples",
    fill = "Amplicon class",
    title = "Amplicon class distribution by therapy timepoint"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

therapy_ecdna_data <- therapy_data %>%
  dplyr::filter(!is.na(patient))

fit_therapy_ecdna <- lme4::glmer(
  ecDNA_present ~ therapy_timepoint + (1 | patient),
  family = binomial,
  data = therapy_ecdna_data
)

beta_therapy <- lme4::fixef(fit_therapy_ecdna)["therapy_timepointPOST-THERAPY"]
se_therapy <- sqrt(
  vcov(fit_therapy_ecdna)[
    "therapy_timepointPOST-THERAPY",
    "therapy_timepointPOST-THERAPY"
  ]
)
p_therapy <- summary(fit_therapy_ecdna)$coefficients[
  "therapy_timepointPOST-THERAPY",
  "Pr(>|z|)"
]

therapy_ecdna_test <- data.frame(
  comparison = "Post-therapy vs Pre-therapy",
  OR = exp(beta_therapy),
  CI_low = exp(beta_therapy - 1.96 * se_therapy),
  CI_high = exp(beta_therapy + 1.96 * se_therapy),
  p_value = p_therapy
)

therapy_ecdna_label <- paste0(
  "OR = ", round(therapy_ecdna_test$OR, 2),
  "\nP = ", signif(therapy_ecdna_test$p_value, 3)
)

therapy_ecdna_plot_data <- therapy_ecdna_data %>%
  dplyr::mutate(
    ecDNA_status = factor(
      ifelse(ecDNA_present == 1, "ecDNA", "non-ecDNA"),
      levels = c("non-ecDNA", "ecDNA")
    )
  )

p_therapy_ecdna <- ggplot2::ggplot(
  therapy_ecdna_plot_data,
  ggplot2::aes(x = therapy_timepoint, fill = ecDNA_status)
) +
  ggplot2::geom_bar(position = "fill", width = 0.60) +
  ggplot2::scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1.10),
    expand = c(0, 0)
  ) +
  ggplot2::annotate(
    "text",
    x = 1.5,
    y = 1.05,
    label = therapy_ecdna_label,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = ecdna_status_colors, drop = FALSE) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = "Therapy timepoint",
    y = "Proportion of samples",
    fill = "Amplicon status",
    title = "ecDNA frequency by therapy timepoint"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold")
  )

# ============================================================
# Recurrent oncogenes by amplicon class, sample, patient, and RMS subtype
# ============================================================

oncogene_parent_class <- cohort_amplicons %>%
  dplyr::mutate(
    Classification = dplyr::case_when(
      ecDNA & FAN ~ "ecDNA + FAN",
      ecDNA ~ "ecDNA",
      FAN ~ "FAN",
      BFB ~ "BFB",
      Complex_non_cyclic ~ "Complex-non-cyclic",
      Linear ~ "Linear",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Classification)) %>%
  dplyr::select(
    patient,
    Sample.name,
    cancer_subtype,
    AA.amplicon.number,
    Classification
  )

oncogene_long <- joined %>%
  dplyr::select(Sample.name, AA.amplicon.number, Oncogenes) %>%
  dplyr::mutate(
    Oncogenes = stringr::str_remove_all(Oncogenes, "\\[|\\]"),
    Oncogenes = stringr::str_remove_all(Oncogenes, "[\"']"),
    Oncogenes = stringr::str_split(Oncogenes, ",\\s*")
  ) %>%
  tidyr::unnest(Oncogenes) %>%
  dplyr::mutate(Oncogenes = stringr::str_trim(Oncogenes)) %>%
  dplyr::filter(!is.na(Oncogenes), Oncogenes != "") %>%
  dplyr::inner_join(
    oncogene_parent_class,
    by = c("Sample.name", "AA.amplicon.number")
  ) %>%
  dplyr::distinct(
    patient,
    Sample.name,
    cancer_subtype,
    AA.amplicon.number,
    Classification,
    Oncogenes
  )

amplicon_levels <- c(
  "ecDNA",
  "ecDNA + FAN",
  "FAN",
  "BFB",
  "Complex-non-cyclic",
  "Linear"
)

amplicon_levels_plot <- amplicon_levels[
  amplicon_levels %in% unique(as.character(oncogene_long$Classification))
]

oncogene_summary_sample <- oncogene_long %>%
  dplyr::distinct(Sample.name, Classification, Oncogenes) %>%
  dplyr::count(Classification, Oncogenes, name = "n")

oncogene_summary_patient <- oncogene_long %>%
  dplyr::distinct(patient, Classification, Oncogenes) %>%
  dplyr::count(Classification, Oncogenes, name = "n")

amplicon_colors <- detailed_amplicon_colors

amplicon_levels <- c(
  "ecDNA",
  "ecDNA + FAN",
  "FAN",
  "BFB",
  "Complex-non-cyclic",
  "Linear"
)

plot_top_oncogenes <- function(summary_df, top_n = 19, y_label = "Count", title = "") {

  top_genes_df <- summary_df %>%
    dplyr::group_by(Oncogenes) %>%
    dplyr::summarise(total_count = sum(n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_count)) %>%
    dplyr::slice_head(n = top_n)

  df_plot <- summary_df %>%
    dplyr::filter(Oncogenes %in% top_genes_df$Oncogenes) %>%
    dplyr::left_join(top_genes_df, by = "Oncogenes") %>%
    dplyr::mutate(
      Oncogenes = reorder(Oncogenes, total_count),
      Classification = factor(Classification, levels = amplicon_levels_plot)
    )

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = Oncogenes, y = n, fill = Classification)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = amplicon_colors,
      drop = TRUE,
      na.translate = FALSE
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::labs(
      title = title,
      x = "Oncogene",
      y = y_label,
      fill = "Amplicon class"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )
}

p_top_oncogene_sample <- plot_top_oncogenes(
  oncogene_summary_sample,
  top_n = 19,
  y_label = "Number of samples"
)

p_top_oncogene_patient <- plot_top_oncogenes(
  oncogene_summary_patient,
  top_n = 19,
  y_label = "Number of patients"
)

# ------------------------------------------------------------
# Recurrent oncogenes by RMS subtype
# ------------------------------------------------------------

plot_oncogenes_by_cancer <- function(
    df,
    level = c("sample", "patient"),
    min_count = 2
) {

  level <- match.arg(level)
  id_col <- ifelse(level == "sample", "Sample.name", "patient")
  y_lab <- ifelse(level == "sample", "Number of samples", "Number of patients")

  subtype_levels <- c("ARMS", "ERMS", "NOS")

  summary_df <- df %>%
    dplyr::filter(cancer_subtype %in% subtype_levels) %>%
    dplyr::distinct(
      .data[[id_col]],
      cancer_subtype,
      Classification,
      Oncogenes
    ) %>%
    dplyr::count(
      cancer_subtype,
      Classification,
      Oncogenes,
      name = "n"
    )

  gene_order <- summary_df %>%
    dplyr::group_by(Oncogenes) %>%
    dplyr::summarise(total = sum(n), .groups = "drop") %>%
    dplyr::filter(total >= min_count) %>%
    dplyr::arrange(total) %>%
    dplyr::pull(Oncogenes)

  plot_df <- summary_df %>%
    dplyr::filter(Oncogenes %in% gene_order) %>%
    tidyr::complete(
      cancer_subtype = subtype_levels,
      Oncogenes = gene_order,
      Classification = amplicon_levels_plot,
      fill = list(n = 0)
    ) %>%
    dplyr::mutate(
      Oncogenes = factor(Oncogenes, levels = gene_order),
      cancer_subtype = factor(cancer_subtype, levels = subtype_levels),
      Classification = factor(Classification, levels = amplicon_levels_plot)
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = n, y = Oncogenes, fill = Classification)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ cancer_subtype) +
    ggplot2::scale_fill_manual(
      values = amplicon_colors,
      drop = TRUE
    ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
    ggplot2::labs(
      x = y_lab,
      y = "Oncogene",
      fill = "Amplicon class"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9, face = "italic")
    )
}
p_cancer_patient <- plot_oncogenes_by_cancer(
  oncogene_long,
  level = "patient",
  min_count = 2
)

p_cancer_sample <- plot_oncogenes_by_cancer(
  oncogene_long,
  level = "sample",
  min_count = 2
)

# ------------------------------------------------------------
# Per-gene view: amplicon class x RMS subtype
# ------------------------------------------------------------

gene_plot_level <- "patient"
gene_id_col <- ifelse(
  gene_plot_level == "sample",
  "Sample.name",
  "patient"
)
gene_y_lab <- ifelse(
  gene_plot_level == "sample",
  "Number of samples",
  "Number of patients"
)

df_gene_plot <- oncogene_long %>%
  dplyr::filter(
    cancer_subtype %in% c("ARMS", "ERMS", "NOS"),
    !is.na(Classification)
  ) %>%
  dplyr::distinct(
    .data[[gene_id_col]],
    Oncogenes,
    Classification,
    cancer_subtype
  ) %>%
  dplyr::count(
    Oncogenes,
    Classification,
    cancer_subtype,
    name = "n"
  )

recurrent_genes <- df_gene_plot %>%
  dplyr::group_by(Oncogenes) %>%
  dplyr::summarise(total = sum(n), .groups = "drop") %>%
  dplyr::filter(total > 1) %>%
  dplyr::arrange(dplyr::desc(total)) %>%
  dplyr::pull(Oncogenes)

df_gene_plot <- df_gene_plot %>%
  dplyr::filter(Oncogenes %in% recurrent_genes) %>%
  dplyr::mutate(
    Oncogenes = factor(Oncogenes, levels = recurrent_genes),
    Classification = factor(
      Classification,
      levels = amplicon_levels
    ),
    cancer_subtype = factor(
      cancer_subtype,
      levels = c("ARMS", "ERMS", "NOS")
    )
  )

subtype_colors <- c(
  "ARMS" = "#2E73B8",
  "ERMS" = "#E69F00",
  "NOS" = "red4"
)

p_gene_subtype <- ggplot2::ggplot(
  df_gene_plot,
  ggplot2::aes(
    x = Classification,
    y = n,
    fill = cancer_subtype
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ Oncogenes) +
  ggplot2::scale_fill_manual(values = subtype_colors) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::labs(
    title = paste0(
      "Distribution of RMS subtypes per oncogene, per ",
      gene_plot_level
    ),
    x = "Amplicon class",
    y = gene_y_lab,
    fill = "RMS subtype"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold", size = 14)
  )

# ============================================================
# ecDNA context: summary and statistics
# ============================================================
# Summarize the ecDNA-context calls already generated above.
ecdna_context_summary <- ecdna_context %>%
  dplyr::count(gene, ecDNA.context, name = "n") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  dplyr::ungroup()

ecdna_context_summary_sample <- ecdna_context %>%
  dplyr::distinct(Sample.name, gene, ecDNA.context) %>%
  dplyr::count(gene, ecDNA.context, name = "n_samples")

ecdna_context_summary_patient <- ecdna_context %>%
  dplyr::filter(!is.na(patient)) %>%
  dplyr::distinct(patient, gene, ecDNA.context) %>%
  dplyr::count(gene, ecDNA.context, name = "n_patients")

ecdna_context_table <- table(
  ecdna_context$gene,
  ecdna_context$ecDNA.context
)

if (
  nrow(ecdna_context_table) >= 2 &&
  ncol(ecdna_context_table) >= 2
) {
  ecdna_context_fisher <- fisher.test(
    ecdna_context_table,
    simulate.p.value = TRUE,
    B = 100000
  )

  ecdna_context_test <- data.frame(
    comparison = "Gene versus ecDNA context",
    test = "Fisher exact test",
    p_value = ecdna_context_fisher$p.value
  )
} else {
  ecdna_context_test <- data.frame(
    comparison = "Gene versus ecDNA context",
    test = "Fisher exact test",
    p_value = NA_real_
  )
}

# ============================================================
# Export key statistical tests together
# ============================================================

key_statistics <- dplyr::bind_rows(
  complexity_test %>%
    dplyr::transmute(
      analysis = "Complexity across amplicon classes",
      statistic = chi_square,
      df = df,
      p_value = p_value
    ),
  oncogene_global_test %>%
    dplyr::transmute(
      analysis = "Oncogene enrichment across amplicon classes",
      statistic = chi_square,
      df = df,
      p_value = p_value
    ),
  oncogene_complexity_test %>%
    dplyr::transmute(
      analysis = "Complexity: oncogene present vs absent",
      statistic = chi_square,
      df = df,
      p_value = p_value
    ),
  copy_number_test %>%
    dplyr::transmute(
      analysis = "Copy number: ecDNA vs chromosomal",
      statistic = chi_square,
      df = df,
      p_value = p_value
    ),
  copy_number_class_global_test %>%
    dplyr::transmute(
      analysis = "Copy number across detailed amplicon classes",
      statistic = chi_square,
      df = df,
      p_value = p_value
    )
)

# ============================================================
# Total ecDNA context counts
# ============================================================

ecdna_context_total <- analysis_features %>%
  dplyr::filter(
    PrimaryFeatureClass == "ecDNA",
    !is.na(ecDNA.context),
    ecDNA.context != "",
    !is.na(Feature.ID)
  ) %>%
  dplyr::distinct(
    patient,
    Sample.name,
    AA.amplicon.number,
    Feature.ID,
    ecDNA.context
  ) %>%
  dplyr::count(ecDNA.context, name = "n_ecDNA_features", sort = TRUE) %>%
  dplyr::mutate(
    total_ecDNA_features_with_context = sum(n_ecDNA_features),
    proportion = n_ecDNA_features / total_ecDNA_features_with_context
  )

p_context_total <- ggplot2::ggplot(
  ecdna_context_total,
  ggplot2::aes(
    x = reorder(ecDNA.context, n_ecDNA_features),
    y = n_ecDNA_features,
    fill = ecDNA.context
  )
) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_ecDNA_features),
    hjust = -0.15,
    size = 4
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggsci::scale_fill_npg() +
  ggplot2::scale_y_continuous(
    breaks = scales::pretty_breaks(),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = "ecDNA context",
    y = "Number of ecDNA features"
  )

# ============================================================
# Pairwise ecDNA-context comparison: MDM2 vs CDK4
# ============================================================

mdm2_cdk4_context <- ecdna_context %>%
  dplyr::filter(gene %in% c("MDM2", "CDK4")) %>%
  dplyr::distinct(
    patient,
    Sample.name,
    Feature.ID,
    gene,
    ecDNA.context
  )

mdm2_cdk4_context_table <- table(
  mdm2_cdk4_context$gene,
  mdm2_cdk4_context$ecDNA.context
)

if (
  nrow(mdm2_cdk4_context_table) == 2 &&
  ncol(mdm2_cdk4_context_table) >= 2
) {
  mdm2_cdk4_fisher <- fisher.test(mdm2_cdk4_context_table)

  mdm2_cdk4_context_test <- data.frame(
    comparison = "MDM2 vs CDK4 ecDNA-context distribution",
    test = "Pairwise Fisher exact test",
    p_value = mdm2_cdk4_fisher$p.value
  )
} else {
  mdm2_cdk4_context_test <- data.frame(
    comparison = "MDM2 vs CDK4 ecDNA-context distribution",
    test = "Pairwise Fisher exact test",
    p_value = NA_real_
  )
}

# ============================================================
# All pairwise ecDNA-context comparisons among recurrent genes
# ============================================================

pairwise_context_genes <- sort(unique(ecdna_context$gene))
pairwise_context_pairs <- combn(pairwise_context_genes, 2, simplify = FALSE)

ecdna_context_pairwise_tests <- lapply(pairwise_context_pairs, function(pair) {

  pair_df <- ecdna_context %>%
    dplyr::filter(gene %in% pair) %>%
    dplyr::distinct(
      patient,
      Sample.name,
      Feature.ID,
      gene,
      ecDNA.context
    )

  pair_table <- table(
    pair_df$gene,
    pair_df$ecDNA.context
  )

  pair_table <- pair_table[
    ,
    colSums(pair_table) > 0,
    drop = FALSE
  ]

  p_value <- if (
    nrow(pair_table) == 2 &&
    ncol(pair_table) >= 2
  ) {
    fisher.test(pair_table)$p.value
  } else {
    NA_real_
  }

  data.frame(
    gene1 = pair[1],
    gene2 = pair[2],
    p_value = p_value
  )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    FDR = p.adjust(p_value, method = "BH"),
    significance = dplyr::case_when(
      FDR < 0.0001 ~ "****",
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  dplyr::arrange(FDR, p_value)

ecdna_context_pairwise_significant <- ecdna_context_pairwise_tests %>%
  dplyr::filter(!is.na(FDR), FDR < 0.05)

# ============================================================
# ecDNA structural context across ALL annotated oncogenes
# ============================================================

ecdna_oncogene_context_all <- analysis_features %>%
  dplyr::filter(
    PrimaryFeatureClass == "ecDNA",
    !is.na(ecDNA.context),
    ecDNA.context != ""
  ) %>%
  dplyr::mutate(
    Oncogenes = clean_list_column(Oncogenes)
  ) %>%
  tidyr::separate_rows(Oncogenes, sep = ",") %>%
  dplyr::mutate(
    oncogene = stringr::str_trim(Oncogenes)
  ) %>%
  dplyr::filter(
    !is.na(oncogene),
    oncogene != "",
    oncogene != "NA",
    oncogene != "NULL"
  ) %>%
  dplyr::distinct(
    patient,
    Sample.name,
    AA.amplicon.number,
    Feature.ID,
    oncogene,
    ecDNA.context
  )

ecdna_oncogene_context_summary_all <- ecdna_oncogene_context_all %>%
  dplyr::count(oncogene, ecDNA.context, name = "n") %>%
  dplyr::group_by(oncogene) %>%
  dplyr::mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  dplyr::ungroup()

ecdna_oncogene_context_table_all <- table(
  ecdna_oncogene_context_all$oncogene,
  ecdna_oncogene_context_all$ecDNA.context
)

if (
  nrow(ecdna_oncogene_context_table_all) >= 2 &&
  ncol(ecdna_oncogene_context_table_all) >= 2
) {
  ecdna_oncogene_context_fisher_all <- fisher.test(
    ecdna_oncogene_context_table_all,
    simulate.p.value = TRUE,
    B = 100000
  )

  ecdna_oncogene_context_global_test_all <- data.frame(
    comparison = "All annotated oncogenes versus ecDNA context",
    test = "Fisher exact test",
    p_value = ecdna_oncogene_context_fisher_all$p.value
  )
} else {
  ecdna_oncogene_context_global_test_all <- data.frame(
    comparison = "All annotated oncogenes versus ecDNA context",
    test = "Fisher exact test",
    p_value = NA_real_
  )
}

recurrent_ecdna_oncogenes_for_context <- ecdna_oncogene_context_all %>%
  dplyr::group_by(oncogene) %>%
  dplyr::summarise(
    n_ecDNA_features = dplyr::n_distinct(Feature.ID),
    n_samples = dplyr::n_distinct(Sample.name),
    n_patients = dplyr::n_distinct(patient),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_ecDNA_features >= 2) %>%
  dplyr::arrange(dplyr::desc(n_ecDNA_features), oncogene)

all_context_pairwise_genes <- recurrent_ecdna_oncogenes_for_context$oncogene

if (length(all_context_pairwise_genes) >= 2) {

  all_context_pairs <- combn(
    all_context_pairwise_genes,
    2,
    simplify = FALSE
  )

  ecdna_oncogene_context_pairwise_all <- lapply(
    all_context_pairs,
    function(pair) {

      pair_df <- ecdna_oncogene_context_all %>%
        dplyr::filter(oncogene %in% pair)

      pair_table <- table(
        pair_df$oncogene,
        pair_df$ecDNA.context
      )

      pair_table <- pair_table[
        ,
        colSums(pair_table) > 0,
        drop = FALSE
      ]

      p_value <- if (
        nrow(pair_table) == 2 &&
        ncol(pair_table) >= 2
      ) {
        fisher.test(pair_table)$p.value
      } else {
        NA_real_
      }

      data.frame(
        oncogene1 = pair[1],
        oncogene2 = pair[2],
        p_value = p_value
      )
    }
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      FDR = p.adjust(p_value, method = "BH"),
      significance = dplyr::case_when(
        FDR < 0.0001 ~ "****",
        FDR < 0.001 ~ "***",
        FDR < 0.01 ~ "**",
        FDR < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::arrange(FDR, p_value)

} else {

  ecdna_oncogene_context_pairwise_all <- data.frame(
    oncogene1 = character(),
    oncogene2 = character(),
    p_value = numeric(),
    FDR = numeric(),
    significance = character()
  )
}

p_all_oncogene_context <- ggplot2::ggplot(
  ecdna_oncogene_context_summary_all,
  ggplot2::aes(
    x = oncogene,
    y = n,
    fill = ecDNA.context
  )
) +
  ggplot2::geom_col(position = "fill") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggsci::scale_fill_npg() +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(
    x = "Oncogene",
    y = "Proportion",
    fill = "ecDNA context"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      face = "italic"
    )
  )

# ============================================================
# Proportion of all annotated oncogene occurrences carried on ecDNA
# ecDNA takes priority for overlapping ecDNA + FAN amplicons
# ============================================================

oncogene_amplicon_long <- joined |>
  dplyr::filter(
    !is.na(Sample.name),
    !is.na(AA.amplicon.number),
    focal_amplicon
  ) |>
  dplyr::mutate(
    Oncogenes = clean_list_column(Oncogenes)
  ) |>
  tidyr::separate_rows(Oncogenes, sep = ",") |>
  dplyr::mutate(
    oncogene = stringr::str_trim(Oncogenes)
  ) |>
  dplyr::filter(
    oncogene != "",
    !is.na(oncogene)
  ) |>
  dplyr::group_by(
    patient,
    Sample.name,
    AA.amplicon.number,
    oncogene
  ) |>
  dplyr::summarise(
    Amplicon_class = dplyr::case_when(
      any(ecDNA_amplicon, na.rm = TRUE) ~ "ecDNA",
      TRUE ~ "Chromosomal"
    ),
    .groups = "drop"
  )

oncogene_ecDNA_proportion <- oncogene_amplicon_long |>
  dplyr::summarise(
    total_oncogene_occurrences = dplyr::n(),
    oncogene_occurrences_on_ecDNA =
      sum(Amplicon_class == "ecDNA"),
    proportion_on_ecDNA =
      oncogene_occurrences_on_ecDNA /
      total_oncogene_occurrences,
    percent_on_ecDNA =
      100 * proportion_on_ecDNA
  )

oncogene_ecDNA_proportion
