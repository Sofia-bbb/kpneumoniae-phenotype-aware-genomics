#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Repository root
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))

input_dir <- file.path(base_dir, "data")
table_dir <- file.path(base_dir, "results", "siderophore")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

side_file <- file.path(input_dir, "04_Siderophore_Table.tsv")
amr_file <- file.path(input_dir, "02_AMR_Table.tsv")
phenotype_file <- file.path(
  input_dir,
  "16_concordance_detail_deduplicated.tsv"
)

side <- read_tsv(
  side_file,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "nan")
)

amr <- read_tsv(
  amr_file,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "nan")
)

phenotype <- read_tsv(
  phenotype_file,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "nan")
)

required_side <- c(
  "SRR",
  "yersiniabactin_present",
  "aerobactin_present",
  "salmochelin_present",
  "siderophore_count",
  "siderophore_profile",
  "virulence_score",
  "resistance_score"
)

required_amr <- c(
  "SRR",
  "num_resistance_classes",
  "num_resistance_genes"
)

missing_side <- setdiff(required_side, names(side))
missing_amr <- setdiff(required_amr, names(amr))

if (length(missing_side) > 0) {
  stop(
    "Missing siderophore columns: ",
    paste(missing_side, collapse = ", ")
  )
}

if (length(missing_amr) > 0) {
  stop(
    "Missing AMR columns: ",
    paste(missing_amr, collapse = ", ")
  )
}

isolate_data <- side |>
  select(all_of(required_side)) |>
  left_join(
    amr |>
      select(all_of(required_amr)),
    by = "SRR"
  ) |>
  mutate(
    across(
      c(
        yersiniabactin_present,
        aerobactin_present,
        salmochelin_present,
        siderophore_count,
        virulence_score,
        resistance_score,
        num_resistance_classes,
        num_resistance_genes
      ),
      as.numeric
    ),

    # Genomic multidrug-resistance burden:
    # determinants affecting at least three resistance classes.
    genomic_MDR = if_else(
      !is.na(num_resistance_classes) &
        num_resistance_classes >= 3,
      1L,
      0L
    )
  )

write_tsv(
  isolate_data,
  file.path(table_dir, "01_isolate_association_dataset.tsv")
)

system_map <- tibble(
  system = c(
    "Yersiniabactin",
    "Aerobactin",
    "Salmochelin"
  ),
  variable = c(
    "yersiniabactin_present",
    "aerobactin_present",
    "salmochelin_present"
  )
)

safe_fisher <- function(a, b, c, d) {
  matrix_data <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE
  )

  test <- fisher.test(matrix_data)

  tibble(
    odds_ratio = unname(test$estimate),
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2],
    p_value = test$p.value
  )
}

# ------------------------------------------------------------
# A. Siderophore systems versus genomic MDR burden
# ------------------------------------------------------------

mdr_results <- lapply(
  seq_len(nrow(system_map)),
  function(i) {

    system_name <- system_map$system[i]
    variable_name <- system_map$variable[i]

    system_present <- isolate_data[[variable_name]] == 1
    mdr_present <- isolate_data$genomic_MDR == 1

    a <- sum(system_present & mdr_present, na.rm = TRUE)
    b <- sum(system_present & !mdr_present, na.rm = TRUE)
    c <- sum(!system_present & mdr_present, na.rm = TRUE)
    d <- sum(!system_present & !mdr_present, na.rm = TRUE)

    stats <- safe_fisher(a, b, c, d)

    tibble(
      siderophore_system = system_name,
      siderophore_MDR = a,
      siderophore_nonMDR = b,
      no_siderophore_MDR = c,
      no_siderophore_nonMDR = d
    ) |>
      bind_cols(stats)
  }
) |>
  bind_rows() |>
  mutate(
    p_adjusted_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    significant_FDR_0_05 = (
      p_adjusted_BH < 0.05
    )
  )

write_tsv(
  mdr_results,
  file.path(
    table_dir,
    "02_siderophore_vs_genomic_MDR.tsv"
  )
)

# ------------------------------------------------------------
# B. Resistance and virulence scores by siderophore system
# ------------------------------------------------------------

continuous_results <- lapply(
  seq_len(nrow(system_map)),
  function(i) {

    system_name <- system_map$system[i]
    variable_name <- system_map$variable[i]

    present <- isolate_data |>
      filter(.data[[variable_name]] == 1)

    absent <- isolate_data |>
      filter(.data[[variable_name]] == 0)

    resistance_test <- wilcox.test(
      present$resistance_score,
      absent$resistance_score,
      exact = FALSE
    )

    virulence_test <- wilcox.test(
      present$virulence_score,
      absent$virulence_score,
      exact = FALSE
    )

    bind_rows(
      tibble(
        siderophore_system = system_name,
        outcome = "resistance_score",
        present_n = sum(
          !is.na(present$resistance_score)
        ),
        absent_n = sum(
          !is.na(absent$resistance_score)
        ),
        present_median = median(
          present$resistance_score,
          na.rm = TRUE
        ),
        absent_median = median(
          absent$resistance_score,
          na.rm = TRUE
        ),
        p_value = resistance_test$p.value
      ),
      tibble(
        siderophore_system = system_name,
        outcome = "virulence_score",
        present_n = sum(
          !is.na(present$virulence_score)
        ),
        absent_n = sum(
          !is.na(absent$virulence_score)
        ),
        present_median = median(
          present$virulence_score,
          na.rm = TRUE
        ),
        absent_median = median(
          absent$virulence_score,
          na.rm = TRUE
        ),
        p_value = virulence_test$p.value
      )
    )
  }
) |>
  bind_rows() |>
  mutate(
    p_adjusted_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    significant_FDR_0_05 = (
      p_adjusted_BH < 0.05
    )
  )

write_tsv(
  continuous_results,
  file.path(
    table_dir,
    "03_siderophore_vs_scores.tsv"
  )
)

# ------------------------------------------------------------
# C. Siderophore systems versus observed phenotype
# ------------------------------------------------------------

phenotype_clean <- phenotype |>
  filter(
    Phenotype %in% c("R", "S")
  ) |>
  select(
    SRR,
    Antibiotic,
    Phenotype
  ) |>
  distinct(
    SRR,
    Antibiotic,
    Phenotype
  ) |>
  left_join(
    isolate_data |>
      select(
        SRR,
        yersiniabactin_present,
        aerobactin_present,
        salmochelin_present,
        siderophore_profile,
        siderophore_count
      ),
    by = "SRR"
  )

phenotype_results <- list()

result_index <- 1

for (antibiotic_name in sort(
  unique(phenotype_clean$Antibiotic)
)) {

  antibiotic_data <- phenotype_clean |>
    filter(Antibiotic == antibiotic_name)

  for (i in seq_len(nrow(system_map))) {

    system_name <- system_map$system[i]
    variable_name <- system_map$variable[i]

    system_present <- (
      antibiotic_data[[variable_name]] == 1
    )

    resistant <- (
      antibiotic_data$Phenotype == "R"
    )

    a <- sum(
      system_present & resistant,
      na.rm = TRUE
    )

    b <- sum(
      system_present & !resistant,
      na.rm = TRUE
    )

    c <- sum(
      !system_present & resistant,
      na.rm = TRUE
    )

    d <- sum(
      !system_present & !resistant,
      na.rm = TRUE
    )

    stats <- safe_fisher(a, b, c, d)

    phenotype_results[[result_index]] <- tibble(
      Antibiotic = antibiotic_name,
      siderophore_system = system_name,
      siderophore_resistant = a,
      siderophore_susceptible = b,
      no_siderophore_resistant = c,
      no_siderophore_susceptible = d,
      evaluated = a + b + c + d
    ) |>
      bind_cols(stats)

    result_index <- result_index + 1
  }
}

phenotype_results <- bind_rows(
  phenotype_results
) |>
  mutate(
    p_adjusted_BH = p.adjust(
      p_value,
      method = "BH"
    ),
    significant_FDR_0_05 = (
      p_adjusted_BH < 0.05
    )
  )

write_tsv(
  phenotype_results,
  file.path(
    table_dir,
    "04_siderophore_vs_antibiotic_phenotype.tsv"
  )
)

# ------------------------------------------------------------
# D. Siderophore profile descriptive integration table
# ------------------------------------------------------------

profile_summary <- isolate_data |>
  group_by(siderophore_profile) |>
  summarise(
    isolates = n(),
    percent = 100 * n() / nrow(isolate_data),
    median_resistance_score = median(
      resistance_score,
      na.rm = TRUE
    ),
    median_virulence_score = median(
      virulence_score,
      na.rm = TRUE
    ),
    genomic_MDR_n = sum(
      genomic_MDR == 1,
      na.rm = TRUE
    ),
    genomic_MDR_percent = 100 * mean(
      genomic_MDR == 1,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(desc(isolates))

write_tsv(
  profile_summary,
  file.path(
    table_dir,
    "05_siderophore_profile_integrated_summary.tsv"
  )
)

cat(
  paste0(
    "\n",
    paste(rep("=", 72), collapse = ""),
    "\n"
  )
)

cat("GPS-Kp Module 9.1: Siderophore Association Analysis\n")

cat(
  paste0(
    paste(rep("=", 72), collapse = ""),
    "\n"
  )
)

cat("Isolates analyzed:", nrow(isolate_data), "\n")
cat("Siderophore systems tested: 3\n")
cat(
  "Antibiotics tested:",
  length(unique(phenotype_clean$Antibiotic)),
  "\n"
)

cat(
  "Significant MDR associations:",
  sum(
    mdr_results$significant_FDR_0_05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant score associations:",
  sum(
    continuous_results$significant_FDR_0_05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Significant phenotype associations:",
  sum(
    phenotype_results$significant_FDR_0_05,
    na.rm = TRUE
  ),
  "\n"
)

cat("Tables written to:", table_dir, "\n")

cat(
  paste0(
    paste(rep("=", 72), collapse = ""),
    "\n"
  )
)

