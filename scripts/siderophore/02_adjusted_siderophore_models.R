#!/usr/bin/env Rscript

.libPaths(c("~/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

cat("=========================================================\n")
cat("GPS-Kp Module 9.2: Adjusted Association Models\n")
cat("=========================================================\n")

dir.create("results/siderophore",
           recursive = TRUE,
           showWarnings = FALSE)

SIDE_INPUT <- "data/04_Siderophore_Table.tsv"
ST_INPUT <- "data/06_ST_Table.tsv"

side_data <- read_tsv(
  SIDE_INPUT,
  show_col_types = FALSE
)

st_data <- read_tsv(
  ST_INPUT,
  show_col_types = FALSE
) %>%
  select(
    SRR,
    K_locus,
    O_locus
  )

dat <- side_data %>%
  left_join(
    st_data,
    by = "SRR"
  )

###############################################################################
## Basic variables
###############################################################################

dat <- dat %>%
  mutate(

    MDR =
      ifelse(resistance_score >= 2,1,0),

    HighVirulence =
      ifelse(virulence_score >= 3,1,0),

    Yersiniabactin =
      ifelse(yersiniabactin_present==1,1,0),

    Aerobactin =
      ifelse(aerobactin_present==1,1,0),

    Salmochelin =
      ifelse(salmochelin_present==1,1,0),

    ST =
      as.factor(ST),

    K_locus =
      as.factor(K_locus)

  )

###############################################################################
## Remove factor variables with only one level
###############################################################################

drop_single_level <- function(df){

  keep <- sapply(df,function(x){

    if(is.factor(x))
      nlevels(droplevels(x)) > 1
    else
      TRUE

  })

  df[,keep]

}

###############################################################################
## Safe logistic regression
###############################################################################

safe_glm <- function(formula,data){

  fit <- tryCatch(

    glm(
      formula,
      family=binomial(),
      data=data
    ),

    error=function(e) NULL

  )

  fit

}

###############################################################################
## Extract odds ratio, confidence interval, and p-value
###############################################################################

extract_result <- function(
  fit,
  exposure,
  outcome,
  model_name
){

  if(is.null(fit)){

    return(
      tibble(
        model = model_name,
        outcome = outcome,
        exposure = exposure,
        n = NA_integer_,
        events = NA_integer_,
        adjusted_odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        converged = FALSE,
        note = "Model fitting failed"
      )
    )

  }

  model_summary <- summary(fit)$coefficients

  possible_rows <- c(
    exposure,
    paste0(exposure,"1")
  )

  coefficient_row <- possible_rows[
    possible_rows %in% rownames(model_summary)
  ]

  if(length(coefficient_row)==0){

    return(
      tibble(
        model = model_name,
        outcome = outcome,
        exposure = exposure,
        n = nobs(fit),
        events = sum(
          model.response(model.frame(fit)) == 1,
          na.rm = TRUE
        ),
        adjusted_odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        converged = isTRUE(fit$converged),
        note = "Exposure coefficient unavailable"
      )
    )

  }

  coefficient_row <- coefficient_row[1]

  estimate <- model_summary[
    coefficient_row,
    "Estimate"
  ]

  standard_error <- model_summary[
    coefficient_row,
    "Std. Error"
  ]

  p_value <- model_summary[
    coefficient_row,
    "Pr(>|z|)"
  ]

  tibble(
    model = model_name,
    outcome = outcome,
    exposure = exposure,
    n = nobs(fit),
    events = sum(
      model.response(model.frame(fit)) == 1,
      na.rm = TRUE
    ),
    adjusted_odds_ratio = exp(estimate),
    ci_low = exp(
      estimate - 1.96 * standard_error
    ),
    ci_high = exp(
      estimate + 1.96 * standard_error
    ),
    p_value = p_value,
    converged = isTRUE(fit$converged),
    note = ""
  )

}

###############################################################################
## Fit one adjusted model safely
###############################################################################

run_adjusted_model <- function(
  data,
  outcome,
  exposure,
  model_name
){

  candidate_variables <- c(
    outcome,
    exposure,
    "ST",
    "Country",
    "Collection_Year",
    "K_locus"
  )

  model_data <- data %>%
    select(
      any_of(candidate_variables)
    ) %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(.data[[exposure]])
    )

  if(nrow(model_data)==0){

    return(
      tibble(
        model = model_name,
        outcome = outcome,
        exposure = exposure,
        n = 0L,
        events = 0L,
        adjusted_odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        converged = FALSE,
        note = "No analyzable records"
      )
    )

  }

  if(
    dplyr::n_distinct(
      model_data[[outcome]],
      na.rm = TRUE
    ) < 2
  ){

    return(
      tibble(
        model = model_name,
        outcome = outcome,
        exposure = exposure,
        n = nrow(model_data),
        events = sum(
          model_data[[outcome]] == 1,
          na.rm = TRUE
        ),
        adjusted_odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        converged = FALSE,
        note = "Outcome has only one level"
      )
    )

  }

  if(
    dplyr::n_distinct(
      model_data[[exposure]],
      na.rm = TRUE
    ) < 2
  ){

    return(
      tibble(
        model = model_name,
        outcome = outcome,
        exposure = exposure,
        n = nrow(model_data),
        events = sum(
          model_data[[outcome]] == 1,
          na.rm = TRUE
        ),
        adjusted_odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        converged = FALSE,
        note = "Exposure has only one level"
      )
    )

  }

  model_data[[exposure]] <- factor(
    model_data[[exposure]]
  )

  adjustment_terms <- character()

  if(
    "ST" %in% names(model_data) &&
    n_distinct(
      model_data$ST,
      na.rm = TRUE
    ) >= 2
  ){

    common_st <- names(
      sort(
        table(model_data$ST),
        decreasing = TRUE
      )
    )[1:min(
      10,
      length(unique(model_data$ST))
    )]

    model_data <- model_data %>%
      mutate(
        ST_group = ifelse(
          ST %in% common_st,
          as.character(ST),
          "Other_ST"
        ),
        ST_group = factor(ST_group)
      )

    adjustment_terms <- c(
      adjustment_terms,
      "ST_group"
    )

  }

  if(
    "Country" %in% names(model_data) &&
    n_distinct(
      model_data$Country,
      na.rm = TRUE
    ) >= 2
  ){

    common_country <- names(
      sort(
        table(model_data$Country),
        decreasing = TRUE
      )
    )[1:min(
      5,
      length(unique(model_data$Country))
    )]

    model_data <- model_data %>%
      mutate(
        Country_group = ifelse(
          Country %in% common_country,
          as.character(Country),
          "Other_country"
        ),
        Country_group = factor(
          Country_group
        )
      )

    adjustment_terms <- c(
      adjustment_terms,
      "Country_group"
    )

  }

  if(
    "Collection_Year" %in% names(model_data)
  ){

    model_data$Collection_Year <-
      suppressWarnings(
        as.numeric(
          model_data$Collection_Year
        )
      )

    if(
      sum(
        !is.na(model_data$Collection_Year)
      ) >= 20 &&
      n_distinct(
        model_data$Collection_Year,
        na.rm = TRUE
      ) >= 2
    ){

      adjustment_terms <- c(
        adjustment_terms,
        "Collection_Year"
      )

    }

  }

  if(
    "K_locus" %in% names(model_data) &&
    n_distinct(
      model_data$K_locus,
      na.rm = TRUE
    ) >= 2
  ){

    common_k <- names(
      sort(
        table(model_data$K_locus),
        decreasing = TRUE
      )
    )[1:min(
      10,
      length(unique(model_data$K_locus))
    )]

    model_data <- model_data %>%
      mutate(
        K_group = ifelse(
          K_locus %in% common_k,
          as.character(K_locus),
          "Other_K"
        ),
        K_group = factor(K_group)
      )

    adjustment_terms <- c(
      adjustment_terms,
      "K_group"
    )

  }

  model_terms <- c(
    exposure,
    adjustment_terms
  )

  variables_needed <- c(
    outcome,
    model_terms
  )

  model_data <- model_data %>%
    filter(
      if_all(
        all_of(variables_needed),
        ~ !is.na(.x)
      )
    ) %>%
    drop_single_level()

  model_terms <- model_terms[
    model_terms %in% names(model_data)
  ]

  if(
    !exposure %in% model_terms
  ){

    return(
      tibble(
        model = model_name,
        outcome = outcome,
        exposure = exposure,
        n = nrow(model_data),
        events = sum(
          model_data[[outcome]] == 1,
          na.rm = TRUE
        ),
        adjusted_odds_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        converged = FALSE,
        note = "Exposure removed after complete-case filtering"
      )
    )

  }

  model_formula <- reformulate(
    termlabels = model_terms,
    response = outcome
  )

  fit <- safe_glm(
    model_formula,
    model_data
  )

  extract_result(
    fit = fit,
    exposure = exposure,
    outcome = outcome,
    model_name = model_name
  )

}

###############################################################################
## Rebuild the analysis dataset using the AMR table
###############################################################################

AMR_INPUT <- "data/02_AMR_Table.tsv"

amr_data <- read_tsv(
  AMR_INPUT,
  show_col_types = FALSE
) %>%
  select(
    SRR,
    num_resistance_classes,
    num_resistance_genes
  )

dat <- dat %>%
  left_join(
    amr_data,
    by = "SRR"
  ) %>%
  mutate(
    num_resistance_classes =
      suppressWarnings(
        as.numeric(num_resistance_classes)
      ),

    num_resistance_genes =
      suppressWarnings(
        as.numeric(num_resistance_genes)
      ),

    # Genomic multidrug-resistance burden:
    # determinants affecting at least three resistance classes.
    genomic_MDR = ifelse(
      !is.na(num_resistance_classes) &
        num_resistance_classes >= 3,
      1,
      0
    )
  )

write_tsv(
  dat,
  "results/siderophore/06_adjusted_model_dataset.tsv"
)

###############################################################################
## Siderophore exposure definitions
###############################################################################

siderophore_map <- tibble(
  siderophore_system = c(
    "Yersiniabactin",
    "Aerobactin",
    "Salmochelin"
  ),

  exposure = c(
    "Yersiniabactin",
    "Aerobactin",
    "Salmochelin"
  )
)

###############################################################################
## A. Adjusted siderophore versus genomic MDR models
###############################################################################

adjusted_mdr_results <- lapply(
  seq_len(nrow(siderophore_map)),
  function(i){

    system_name <-
      siderophore_map$siderophore_system[i]

    exposure_name <-
      siderophore_map$exposure[i]

    run_adjusted_model(
      data = dat,
      outcome = "genomic_MDR",
      exposure = exposure_name,
      model_name = "Adjusted genomic MDR model"
    ) %>%
      mutate(
        siderophore_system = system_name
      )

  }
) %>%
  bind_rows() %>%
  mutate(
    p_adjusted_BH = p.adjust(
      p_value,
      method = "BH"
    ),

    significant_FDR_0_05 =
      !is.na(p_adjusted_BH) &
      p_adjusted_BH < 0.05
  ) %>%
  select(
    siderophore_system,
    everything()
  )

write_tsv(
  adjusted_mdr_results,
  "results/siderophore/07_adjusted_siderophore_vs_genomic_MDR.tsv"
)

###############################################################################
## B. Prepare observed antibiotic phenotype data
###############################################################################

PHENOTYPE_INPUT <-
  "data/16_concordance_detail_deduplicated.tsv"

phenotype_raw <- read_tsv(
  PHENOTYPE_INPUT,
  show_col_types = FALSE
)

phenotype_data <- phenotype_raw %>%
  filter(
    Phenotype %in% c("R", "S")
  ) %>%
  transmute(
    SRR,
    Antibiotic,

    resistant = ifelse(
      Phenotype == "R",
      1,
      0
    )
  ) %>%
  distinct(
    SRR,
    Antibiotic,
    .keep_all = TRUE
  ) %>%
  left_join(
    dat %>%
      select(
        SRR,
        Yersiniabactin,
        Aerobactin,
        Salmochelin,
        ST,
        Country,
        Collection_Year,
        K_locus
      ),
    by = "SRR"
  )

write_tsv(
  phenotype_data,
  "results/siderophore/08_adjusted_phenotype_model_dataset.tsv"
)

###############################################################################
## C. Adjusted siderophore versus antibiotic phenotype models
###############################################################################

phenotype_results <- list()
result_number <- 1

antibiotic_names <- sort(
  unique(phenotype_data$Antibiotic)
)

for(antibiotic_name in antibiotic_names){

  antibiotic_data <- phenotype_data %>%
    filter(
      Antibiotic == antibiotic_name
    )

  for(i in seq_len(nrow(siderophore_map))){

    system_name <-
      siderophore_map$siderophore_system[i]

    exposure_name <-
      siderophore_map$exposure[i]

    current_result <- run_adjusted_model(
      data = antibiotic_data,
      outcome = "resistant",
      exposure = exposure_name,
      model_name = "Adjusted phenotype resistance model"
    ) %>%
      mutate(
        Antibiotic = antibiotic_name,
        siderophore_system = system_name
      )

    phenotype_results[[result_number]] <-
      current_result

    result_number <- result_number + 1
  }
}

adjusted_phenotype_results <- bind_rows(
  phenotype_results
) %>%
  mutate(
    p_adjusted_BH = p.adjust(
      p_value,
      method = "BH"
    ),

    significant_FDR_0_05 =
      !is.na(p_adjusted_BH) &
      p_adjusted_BH < 0.05
  ) %>%
  select(
    Antibiotic,
    siderophore_system,
    everything()
  )

write_tsv(
  adjusted_phenotype_results,
  "results/siderophore/08_adjusted_siderophore_vs_phenotype.tsv"
)

###############################################################################
## D. Compare adjusted and unadjusted genomic MDR estimates
###############################################################################

UNADJUSTED_MDR <-
  "results/siderophore/02_siderophore_vs_genomic_MDR.tsv"

if(file.exists(UNADJUSTED_MDR)){

  unadjusted_mdr <- read_tsv(
    UNADJUSTED_MDR,
    show_col_types = FALSE
  ) %>%
    select(
      siderophore_system,

      unadjusted_odds_ratio =
        odds_ratio,

      unadjusted_ci_low =
        ci_low,

      unadjusted_ci_high =
        ci_high,

      unadjusted_p_adjusted_BH =
        p_adjusted_BH
    )

  mdr_comparison <- adjusted_mdr_results %>%
    left_join(
      unadjusted_mdr,
      by = "siderophore_system"
    )

  write_tsv(
    mdr_comparison,
    "results/siderophore/09_MDR_unadjusted_vs_adjusted.tsv"
  )
}

###############################################################################
## E. Compare adjusted and unadjusted phenotype estimates
###############################################################################

UNADJUSTED_PHENOTYPE <-
  "results/siderophore/04_siderophore_vs_antibiotic_phenotype.tsv"

if(file.exists(UNADJUSTED_PHENOTYPE)){

  unadjusted_phenotype <- read_tsv(
    UNADJUSTED_PHENOTYPE,
    show_col_types = FALSE
  ) %>%
    select(
      Antibiotic,
      siderophore_system,

      unadjusted_odds_ratio =
        odds_ratio,

      unadjusted_ci_low =
        ci_low,

      unadjusted_ci_high =
        ci_high,

      unadjusted_p_adjusted_BH =
        p_adjusted_BH
    )

  phenotype_comparison <-
    adjusted_phenotype_results %>%
    left_join(
      unadjusted_phenotype,
      by = c(
        "Antibiotic",
        "siderophore_system"
      )
    )

  write_tsv(
    phenotype_comparison,
    "results/siderophore/10_phenotype_unadjusted_vs_adjusted.tsv"
  )
}

###############################################################################
## Final summary
###############################################################################

cat("\n")
cat("=========================================================\n")
cat("GPS-Kp Module 9.2 completed\n")
cat("=========================================================\n")

cat(
  "Isolates analyzed:",
  nrow(dat),
  "\n"
)

cat(
  "Adjusted genomic MDR models:",
  nrow(adjusted_mdr_results),
  "\n"
)

cat(
  "Significant adjusted MDR associations:",
  sum(
    adjusted_mdr_results$significant_FDR_0_05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Adjusted phenotype models:",
  nrow(adjusted_phenotype_results),
  "\n"
)

cat(
  "Significant adjusted phenotype associations:",
  sum(
    adjusted_phenotype_results$significant_FDR_0_05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Failed or non-converged MDR models:",
  sum(
    !adjusted_mdr_results$converged,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Failed or non-converged phenotype models:",
  sum(
    !adjusted_phenotype_results$converged,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Results directory:",
  "results/siderophore",
  "\n"
)

cat("=========================================================\n")



