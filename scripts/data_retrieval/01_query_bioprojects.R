#!/usr/bin/env Rscript

# ============================================================
# Galaxy-NUBRI NCBI Data Retrieval Pipeline
# Module: 01_query_bioprojects.R
#
# Purpose:
#   Read BioProject accessions extracted from the PubMed
#   pipeline and count linked BioSample and SRA records.
#
# Input:
#   outputs/04_accessions_long.csv
#
# Outputs:
#   NCBI_pipeline/outputs/01_bioproject_summary.csv
#   NCBI_pipeline/logs/01_query_bioprojects_log.txt
#
# Version: 1.0.0
# ============================================================

suppressPackageStartupMessages({
  library(rentrez)
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

input_file <- "outputs/04_accessions_long.csv"

output_file <- paste0(
  "NCBI_pipeline/outputs/",
  "01_bioproject_summary.csv"
)

log_file <- paste0(
  "NCBI_pipeline/logs/",
  "01_query_bioprojects_log.txt"
)

dir.create(
  "NCBI_pipeline/outputs",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "NCBI_pipeline/logs",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# NCBI identification
# ------------------------------------------------------------


# Add an NCBI API key here later if you create one:
# set_entrez_key("YOUR_NCBI_API_KEY")

# ------------------------------------------------------------
# Logging
# ------------------------------------------------------------

log_message <- function(message) {

  line <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    message
  )

  cat(line, "\n")
  cat(
    line,
    "\n",
    file = log_file,
    append = TRUE
  )
}

writeLines(
  c(
    "Galaxy-NUBRI NCBI Data Retrieval Pipeline",
    "Module: 01_query_bioprojects.R",
    "Version: 1.0.0",
    paste("Start time:", Sys.time()),
    ""
  ),
  log_file
)

# ------------------------------------------------------------
# Validate input
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

accession_data <- read_csv(
  input_file,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

required_columns <- c(
  "PMID",
  "PMCID",
  "Accession_Type",
  "Accession"
)

missing_columns <- setdiff(
  required_columns,
  names(accession_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Extract unique BioProjects
# ------------------------------------------------------------

bioprojects <- accession_data |>
  filter(Accession_Type == "BioProject") |>
  mutate(
    BioProject = str_to_upper(str_trim(Accession))
  ) |>
  filter(
    str_detect(
      BioProject,
      "^PRJ(NA|EB|DB)[0-9]+$"
    )
  ) |>
  group_by(BioProject) |>
  summarise(
    PMID = paste(
      sort(unique(PMID)),
      collapse = ";"
    ),
    PMCID = paste(
      sort(unique(PMCID)),
      collapse = ";"
    ),
    Publication_Count = n_distinct(PMID),
    .groups = "drop"
  ) |>
  arrange(BioProject)

total_projects <- nrow(bioprojects)

log_message(
  paste(
    "Unique BioProjects loaded:",
    total_projects
  )
)

if (total_projects == 0) {
  stop("No BioProject accessions were found.")
}

# ------------------------------------------------------------
# Safe Entrez count
# ------------------------------------------------------------

safe_entrez_count <- function(
  database,
  search_term
) {

  result <- tryCatch(
    {
      search_result <- entrez_search(
        db = database,
        term = search_term,
        retmax = 0
      )

      list(
        count = as.integer(search_result$count),
        status = "Success",
        error = NA_character_
      )
    },
    error = function(error_object) {

      list(
        count = NA_integer_,
        status = "Failed",
        error = conditionMessage(error_object)
      )
    }
  )

  result
}

# ------------------------------------------------------------
# Query one BioProject
# ------------------------------------------------------------

query_one_project <- function(
  bioproject,
  pmid,
  pmcid,
  publication_count,
  row_number,
  total_rows
) {

  log_message(
    paste(
      "Querying",
      row_number,
      "of",
      total_rows,
      "-",
      bioproject
    )
  )

  bioproject_query <- paste0(
    bioproject,
    "[BioProject]"
  )

  project_result <- safe_entrez_count(
    database = "bioproject",
    search_term = bioproject_query
  )

  Sys.sleep(0.4)

  biosample_result <- safe_entrez_count(
    database = "biosample",
    search_term = bioproject_query
  )

  Sys.sleep(0.4)

  sra_result <- safe_entrez_count(
    database = "sra",
    search_term = bioproject_query
  )

  Sys.sleep(0.4)

  overall_status <- if (
    project_result$status == "Success" &&
    biosample_result$status == "Success" &&
    sra_result$status == "Success"
  ) {
    "Completed"
  } else {
    "Partial or failed"
  }

  combined_errors <- paste(
    na.omit(
      c(
        project_result$error,
        biosample_result$error,
        sra_result$error
      )
    ),
    collapse = " | "
  )

  if (combined_errors == "") {
    combined_errors <- NA_character_
  }

  tibble(
    BioProject = bioproject,
    PMID = pmid,
    PMCID = pmcid,
    Publication_Count = publication_count,

    BioProject_Record_Count =
      project_result$count,

    BioSample_Record_Count =
      biosample_result$count,

    SRA_Record_Count =
      sra_result$count,

    BioProject_Status =
      project_result$status,

    BioSample_Status =
      biosample_result$status,

    SRA_Status =
      sra_result$status,

    Query_Status =
      overall_status,

    Query_Error =
      combined_errors
  )
}

# ------------------------------------------------------------
# Run all queries
# ------------------------------------------------------------

results <- map_dfr(
  seq_len(total_projects),
  function(i) {

    query_one_project(
      bioproject =
        bioprojects$BioProject[i],

      pmid =
        bioprojects$PMID[i],

      pmcid =
        bioprojects$PMCID[i],

      publication_count =
        bioprojects$Publication_Count[i],

      row_number = i,
      total_rows = total_projects
    )
  }
)

# ------------------------------------------------------------
# Add useful summary fields
# ------------------------------------------------------------

results <- results |>
  mutate(
    Has_BioSamples = case_when(
      is.na(BioSample_Record_Count) ~ NA_character_,
      BioSample_Record_Count > 0 ~ "Yes",
      TRUE ~ "No"
    ),

    Has_SRA_Data = case_when(
      is.na(SRA_Record_Count) ~ NA_character_,
      SRA_Record_Count > 0 ~ "Yes",
      TRUE ~ "No"
    ),

    Retrieval_Priority = case_when(
      BioSample_Record_Count > 0 &
        SRA_Record_Count > 0 ~ "Ready",

      BioSample_Record_Count > 0 &
        SRA_Record_Count == 0 ~
        "BioSamples only",

      BioSample_Record_Count == 0 &
        SRA_Record_Count > 0 ~
        "SRA only",

      BioSample_Record_Count == 0 &
        SRA_Record_Count == 0 ~
        "No linked records",

      TRUE ~ "Check query"
    )
  ) |>
  arrange(
    desc(SRA_Record_Count),
    desc(BioSample_Record_Count),
    BioProject
  )

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------

write_csv(
  results,
  output_file,
  na = ""
)

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

completed_count <- sum(
  results$Query_Status == "Completed",
  na.rm = TRUE
)

ready_count <- sum(
  results$Retrieval_Priority == "Ready",
  na.rm = TRUE
)

total_biosamples <- sum(
  results$BioSample_Record_Count,
  na.rm = TRUE
)

total_sra_records <- sum(
  results$SRA_Record_Count,
  na.rm = TRUE
)

log_message(
  paste(
    "BioProjects processed:",
    nrow(results)
  )
)

log_message(
  paste(
    "Queries completed:",
    completed_count
  )
)

log_message(
  paste(
    "Projects ready for retrieval:",
    ready_count
  )
)

log_message(
  paste(
    "Total linked BioSample records:",
    total_biosamples
  )
)

log_message(
  paste(
    "Total linked SRA records:",
    total_sra_records
  )
)

log_message(
  paste(
    "Output file:",
    output_file
  )
)

log_message("Module completed successfully.")

cat("\n")
cat("NCBI Module 1 completed successfully.\n")
cat(
  "BioProjects processed:",
  nrow(results),
  "\n"
)
cat(
  "Queries completed:",
  completed_count,
  "\n"
)
cat(
  "Projects ready for retrieval:",
  ready_count,
  "\n"
)
cat(
  "Total linked BioSample records:",
  total_biosamples,
  "\n"
)
cat(
  "Total linked SRA records:",
  total_sra_records,
  "\n"
)
cat(
  "Output:",
  output_file,
  "\n"
)

