#!/usr/bin/env Rscript

# ==============================================================================
# Module 4A: AST/MIC Extraction from NCBI BioSample Metadata
# Project: Klebsiella pneumoniae NCBI Data Retrieval Pipeline
#
# Inputs:
#   outputs/03_clean_metadata.csv
#   outputs/03_biosample_metadata.xml
#
# Outputs:
#   outputs/04_ast_mic_candidates.csv
#   outputs/04_ast_attribute_long.csv
#   logs/04_ast_mic_summary.csv
#
# Purpose:
#   Extract antimicrobial susceptibility information recorded in BioSample
#   attributes for:
#       - Meropenem
#       - Imipenem
#       - Ciprofloxacin
#       - Colistin
#
# Important:
#   This script does not infer missing susceptibility results.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(xml2)
})


# ------------------------------------------------------------------------------
# 2. Define paths
# ------------------------------------------------------------------------------

metadata_file <- "outputs/03_clean_metadata.csv"
biosample_xml_file <- "outputs/03_biosample_metadata.xml"

output_file <- "outputs/04_ast_mic_candidates.csv"
long_output_file <- "outputs/04_ast_attribute_long.csv"
summary_file <- "logs/04_ast_mic_summary.csv"


# ------------------------------------------------------------------------------
# 3. Check input files
# ------------------------------------------------------------------------------

if (!file.exists(metadata_file)) {
  stop("Missing input file: ", metadata_file)
}

if (!file.exists(biosample_xml_file)) {
  stop("Missing BioSample XML file: ", biosample_xml_file)
}

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------------------------

clean_missing <- function(x) {

  x <- as.character(x)
  x <- str_squish(x)

  missing_terms <- c(
    "",
    "NA",
    "N/A",
    "na",
    "n/a",
    "null",
    "NULL",
    "unknown",
    "Unknown",
    "not available",
    "Not available",
    "not provided",
    "Not provided",
    "not collected",
    "Not collected",
    "-"
  )

  x[x %in% missing_terms] <- NA_character_

  x
}


normalize_attribute_name <- function(x) {

  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}


identify_antibiotic <- function(attribute_name, attribute_value) {

  combined_text <- paste(
    attribute_name,
    attribute_value,
    sep = " "
  ) %>%
    str_to_lower()

  case_when(
    str_detect(combined_text, "meropenem|\\bmem\\b") ~ "Meropenem",
    str_detect(combined_text, "imipenem|\\bipm\\b") ~ "Imipenem",
    str_detect(combined_text, "ciprofloxacin|\\bcip\\b") ~ "Ciprofloxacin",
    str_detect(combined_text, "colistin|polymyxin[ _-]?e") ~ "Colistin",
    TRUE ~ NA_character_
  )
}


extract_interpretation <- function(x) {

  x <- clean_missing(x)
  x_lower <- str_to_lower(x)

  case_when(
    is.na(x_lower) ~ NA_character_,

    str_detect(
      x_lower,
      "(^|[^a-z])(resistant|resistance|non[- ]?susceptible|\\br\\b)([^a-z]|$)"
    ) ~ "R",

    str_detect(
      x_lower,
      "(^|[^a-z])(intermediate|\\bi\\b)([^a-z]|$)"
    ) ~ "I",

    str_detect(
      x_lower,
      "(^|[^a-z])(susceptible|sensitive|\\bs\\b)([^a-z]|$)"
    ) ~ "S",

    TRUE ~ NA_character_
  )
}


extract_mic <- function(x) {

  x <- clean_missing(x)

  mic_text <- str_extract(
    x,
    regex(
      "(MIC\\s*[:=]?\\s*)?[<>≤≥]?[[:space:]]*[0-9]+(?:\\.[0-9]+)?[[:space:]]*(µg/mL|ug/mL|mg/L)",
      ignore_case = TRUE
    )
  )

  clean_missing(mic_text)
}


combine_result <- function(value, interpretation, mic) {

  pieces <- c(
    clean_missing(value),
    clean_missing(interpretation),
    clean_missing(mic)
  )

  pieces <- unique(pieces[!is.na(pieces)])

  if (length(pieces) == 0) {
    return(NA_character_)
  }

  paste(pieces, collapse = " | ")
}


collapse_unique_values <- function(x) {

  x <- clean_missing(x)
  x <- unique(x[!is.na(x)])

  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = " ; ")
}


# ------------------------------------------------------------------------------
# 5. Read Module 3 metadata
# ------------------------------------------------------------------------------

message("Reading Module 3 metadata...")

metadata <- read_csv(
  metadata_file,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "Run_Accession",
  "BioProject",
  "BioSample",
  "Sample_Name",
  "Country",
  "Collection_Year",
  "Isolation_Source"
)

missing_columns <- setdiff(required_columns, names(metadata))

if (length(missing_columns) > 0) {
  stop(
    "Missing columns in Module 3 metadata: ",
    paste(missing_columns, collapse = ", ")
  )
}

metadata <- metadata %>%
  mutate(across(everything(), clean_missing))


# ------------------------------------------------------------------------------
# 6. Read BioSample XML
# ------------------------------------------------------------------------------

message("Reading BioSample XML...")

biosample_document <- read_xml(biosample_xml_file)

sample_nodes <- xml_find_all(
  biosample_document,
  ".//BioSample"
)

message("BioSample records found: ", length(sample_nodes))


# ------------------------------------------------------------------------------
# 7. Extract all BioSample attributes into long format
# ------------------------------------------------------------------------------

message("Extracting BioSample attributes...")

attribute_long <- map_dfr(
  sample_nodes,
  function(sample_node) {

    biosample_accession <- xml_attr(
      sample_node,
      "accession"
    )

    attribute_nodes <- xml_find_all(
      sample_node,
      ".//Attribute"
    )

    if (length(attribute_nodes) == 0) {
      return(tibble())
    }

    attribute_name <- xml_attr(
      attribute_nodes,
      "harmonized_name"
    )

    original_name <- xml_attr(
      attribute_nodes,
      "attribute_name"
    )

    attribute_name <- ifelse(
      !is.na(attribute_name) & attribute_name != "",
      attribute_name,
      original_name
    )

    tibble(
      BioSample = biosample_accession,
      Attribute_Name = attribute_name,
      Attribute_Value = xml_text(
        attribute_nodes,
        trim = TRUE
      )
    )
  }
) %>%
  mutate(
    BioSample = clean_missing(BioSample),

    Attribute_Name = clean_missing(
      Attribute_Name
    ),

    Attribute_Name_Normalized = normalize_attribute_name(
      Attribute_Name
    ),

    Attribute_Value = clean_missing(
      Attribute_Value
    )
  )


# ------------------------------------------------------------------------------
# 8. Identify AST-related attributes
# ------------------------------------------------------------------------------

ast_keywords <- paste(
  c(
    "antibiotic",
    "antimicrobial",
    "suscept",
    "resistan",
    "mic",
    "minimum_inhibitory",
    "meropenem",
    "imipenem",
    "ciprofloxacin",
    "colistin",
    "polymyxin"
  ),
  collapse = "|"
)

ast_attribute_long <- attribute_long %>%
  filter(
    str_detect(
      str_to_lower(
        paste(
          Attribute_Name_Normalized,
          Attribute_Value
        )
      ),
      ast_keywords
    )
  ) %>%
  mutate(
    Antibiotic = map2_chr(
      Attribute_Name_Normalized,
      Attribute_Value,
      identify_antibiotic
    ),

    Interpretation = map_chr(
      Attribute_Value,
      extract_interpretation
    ),

    MIC = map_chr(
      Attribute_Value,
      extract_mic
    )
  )

write_csv(
  ast_attribute_long,
  long_output_file,
  na = ""
)


# ------------------------------------------------------------------------------
# 9. Retain the four target antibiotics
# ------------------------------------------------------------------------------

target_ast <- ast_attribute_long %>%
  filter(!is.na(Antibiotic)) %>%
  rowwise() %>%
  mutate(
    Standardized_Result = combine_result(
      Attribute_Value,
      Interpretation,
      MIC
    )
  ) %>%
  ungroup()


# ------------------------------------------------------------------------------
# 10. Convert antibiotic results to wide format
# ------------------------------------------------------------------------------

ast_wide <- target_ast %>%
  group_by(
    BioSample,
    Antibiotic
  ) %>%
  summarise(
    Result = collapse_unique_values(
      Standardized_Result
    ),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Antibiotic,
    values_from = Result
  )


expected_antibiotics <- c(
  "Meropenem",
  "Imipenem",
  "Ciprofloxacin",
  "Colistin"
)

for (antibiotic in expected_antibiotics) {

  if (!antibiotic %in% names(ast_wide)) {
    ast_wide[[antibiotic]] <- NA_character_
  }
}


# ------------------------------------------------------------------------------
# 11. Join AST information to Module 3 metadata
# ------------------------------------------------------------------------------

ast_candidates <- metadata %>%
  left_join(
    ast_wide,
    by = "BioSample"
  ) %>%
  mutate(
    Resistance_Phenotype = case_when(
      str_detect(
        paste(
          Meropenem,
          Imipenem,
          Ciprofloxacin,
          Colistin
        ),
        "(^|[^A-Z])R([^A-Z]|$)"
      ) ~ "Resistance reported",

      if_any(
        all_of(expected_antibiotics),
        ~ !is.na(.x)
      ) ~ "AST reported; manual review required",

      TRUE ~ NA_character_
    ),

    AST_Data_Source = case_when(
      if_any(
        all_of(expected_antibiotics),
        ~ !is.na(.x)
      ) ~ "NCBI BioSample",

      TRUE ~ NA_character_
    ),

    Manual_Review_Status = case_when(
      if_any(
        all_of(expected_antibiotics),
        ~ !is.na(.x)
      ) ~ "Review extracted value",

      TRUE ~ "AST not found in BioSample"
    )
  ) %>%
  select(
    Run_Accession,
    BioProject,
    BioSample,
    Sample_Name,
    Country,
    Collection_Year,
    Isolation_Source,
    Meropenem,
    Imipenem,
    Ciprofloxacin,
    Colistin,
    Resistance_Phenotype,
    AST_Data_Source,
    Manual_Review_Status
  )


# ------------------------------------------------------------------------------
# 12. Export candidate AST table
# ------------------------------------------------------------------------------

write_csv(
  ast_candidates,
  output_file,
  na = ""
)


# ------------------------------------------------------------------------------
# 13. Create summary
# ------------------------------------------------------------------------------

ast_present <- rowSums(
  !is.na(
    ast_candidates[, expected_antibiotics]
  )
) > 0

summary_table <- tibble(
  Metric = c(
    "Metadata records imported",
    "Unique BioSamples in metadata",
    "BioSample XML records",
    "Total BioSample attributes extracted",
    "AST-related attributes detected",
    "Records with meropenem information",
    "Records with imipenem information",
    "Records with ciprofloxacin information",
    "Records with colistin information",
    "Records with at least one target antibiotic",
    "Records without target AST information"
  ),

  Value = c(
    nrow(metadata),
    n_distinct(metadata$BioSample),
    length(sample_nodes),
    nrow(attribute_long),
    nrow(ast_attribute_long),
    sum(!is.na(ast_candidates$Meropenem)),
    sum(!is.na(ast_candidates$Imipenem)),
    sum(!is.na(ast_candidates$Ciprofloxacin)),
    sum(!is.na(ast_candidates$Colistin)),
    sum(ast_present),
    sum(!ast_present)
  )
)

write_csv(
  summary_table,
  summary_file
)


# ------------------------------------------------------------------------------
# 14. Completion report
# ------------------------------------------------------------------------------

message("")
message("============================================================")
message("Module 4A completed")
message("============================================================")
message("Metadata records: ", nrow(metadata))
message("AST-related attributes: ", nrow(ast_attribute_long))
message("")
message("Candidate AST table:")
message("  ", output_file)
message("")
message("Long-format AST attributes:")
message("  ", long_output_file)
message("")
message("Summary:")
message("  ", summary_file)
message("============================================================")

