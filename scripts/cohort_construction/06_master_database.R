#!/usr/bin/env Rscript

# ============================================================
# Module 06: Master Genotype-Phenotype Database
# NCBI Pathogen Detection Genotype-Phenotype Pipeline
#
# Purpose:
#   Merge the human clinical sequencing/metadata records with
#   the harmonized AST/MIC data using BioSample accession.
#
# Inputs:
#   outputs/05_Human_Clinical.csv
#   outputs/04_AST_MIC_long.csv
#
# Outputs:
#   outputs/06_Master_Genotype_Phenotype_long.csv
#   outputs/06_Human_Clinical_Run_List.csv
#   outputs/06_Unmatched_Human_Metadata.csv
#   outputs/06_Unmatched_AST_BioSamples.csv
#   logs/06_master_database_summary.csv
#
# Important:
#   The master database remains in LONG format:
#   one row per sequencing run and antimicrobial observation.
# ============================================================


# -----------------------------
# 1. File locations
# -----------------------------

metadata_file <- "outputs/05_Human_Clinical.csv"
ast_file <- "outputs/04_AST_MIC_long.csv"

master_output <- "outputs/06_Master_Genotype_Phenotype_long.csv"
run_list_output <- "outputs/06_Human_Clinical_Run_List.csv"

unmatched_metadata_output <-
  "outputs/06_Unmatched_Human_Metadata.csv"

unmatched_ast_output <-
  "outputs/06_Unmatched_AST_BioSamples.csv"

summary_output <- "logs/06_master_database_summary.csv"


# -----------------------------
# 2. Check input files
# -----------------------------

required_files <- c(
  metadata_file,
  ast_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    paste(
      "Required input file(s) not found:",
      paste(missing_files, collapse = ", ")
    )
  )
}


# -----------------------------
# 3. Read input datasets
# -----------------------------

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c(
    "",
    "NA",
    "N/A",
    "na",
    "n/a",
    "NULL",
    "null"
  )
)

ast <- read.csv(
  ast_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c(
    "",
    "NA",
    "N/A",
    "na",
    "n/a",
    "NULL",
    "null"
  )
)

cat("Human-clinical metadata rows:", nrow(metadata), "\n")
cat("AST/MIC rows:", nrow(ast), "\n")


# -----------------------------
# 4. Normalize names for matching
# -----------------------------

normalize_name <- function(x) {

  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)

  return(x)
}


find_column <- function(data, possible_names) {

  normalized_actual <- normalize_name(names(data))
  normalized_possible <- normalize_name(possible_names)

  match_index <- match(
    normalized_possible,
    normalized_actual,
    nomatch = 0
  )

  match_index <- match_index[
    match_index > 0
  ]

  if (length(match_index) == 0) {
    return(NULL)
  }

  names(data)[match_index[1]]
}


# -----------------------------
# 5. Identify key columns
# -----------------------------

metadata_biosample_column <- find_column(
  metadata,
  c(
    "BioSample",
    "biosample",
    "biosample_accession",
    "BioSample_Accession"
  )
)

ast_biosample_column <- find_column(
  ast,
  c(
    "BioSample",
    "biosample",
    "biosample_accession",
    "BioSample_Accession"
  )
)

run_column <- find_column(
  metadata,
  c(
    "Run_Accession",
    "Run",
    "run_accession",
    "SRA_Run",
    "sra_accession"
  )
)


if (is.null(metadata_biosample_column)) {

  stop(
    "BioSample column was not found in the human-clinical metadata."
  )
}

if (is.null(ast_biosample_column)) {

  stop(
    "BioSample column was not found in the AST/MIC dataset."
  )
}


cat(
  "Metadata BioSample column:",
  metadata_biosample_column,
  "\n"
)

cat(
  "AST BioSample column:",
  ast_biosample_column,
  "\n"
)

cat(
  "Run column:",
  ifelse(
    is.null(run_column),
    "Not found",
    run_column
  ),
  "\n"
)


# -----------------------------
# 6. Standardize BioSample values
# -----------------------------

clean_accession <- function(x) {

  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)

  x[
    x %in% c(
      "",
      "NA",
      "N/A",
      "NULL"
    )
  ] <- NA_character_

  return(x)
}


metadata$BioSample_Merge_Key <- clean_accession(
  metadata[[metadata_biosample_column]]
)

ast$BioSample_Merge_Key <- clean_accession(
  ast[[ast_biosample_column]]
)


# Remove rows without BioSample identifiers

metadata_missing_biosample <- sum(
  is.na(metadata$BioSample_Merge_Key)
)

ast_missing_biosample <- sum(
  is.na(ast$BioSample_Merge_Key)
)

metadata <- metadata[
  !is.na(metadata$BioSample_Merge_Key),
  ,
  drop = FALSE
]

ast <- ast[
  !is.na(ast$BioSample_Merge_Key),
  ,
  drop = FALSE
]


# -----------------------------
# 7. Remove exact duplicate rows
# -----------------------------

metadata_before_dedup <- nrow(metadata)
ast_before_dedup <- nrow(ast)

metadata <- unique(metadata)
ast <- unique(ast)

metadata_duplicates_removed <-
  metadata_before_dedup - nrow(metadata)

ast_duplicates_removed <-
  ast_before_dedup - nrow(ast)


# -----------------------------
# 8. Determine matched BioSamples
# -----------------------------

metadata_biosamples <- unique(
  metadata$BioSample_Merge_Key
)

ast_biosamples <- unique(
  ast$BioSample_Merge_Key
)

matched_biosamples <- intersect(
  metadata_biosamples,
  ast_biosamples
)

metadata_only_biosamples <- setdiff(
  metadata_biosamples,
  ast_biosamples
)

ast_only_biosamples <- setdiff(
  ast_biosamples,
  metadata_biosamples
)


# -----------------------------
# 9. Create unmatched reports
# -----------------------------

unmatched_metadata <- metadata[
  metadata$BioSample_Merge_Key %in%
    metadata_only_biosamples,
  ,
  drop = FALSE
]

unmatched_ast <- ast[
  ast$BioSample_Merge_Key %in%
    ast_only_biosamples,
  ,
  drop = FALSE
]

write.csv(
  unmatched_metadata,
  unmatched_metadata_output,
  row.names = FALSE,
  na = ""
)

write.csv(
  unmatched_ast,
  unmatched_ast_output,
  row.names = FALSE,
  na = ""
)


# -----------------------------
# 10. Merge metadata with AST/MIC
# -----------------------------

master <- merge(
  metadata,
  ast,
  by = "BioSample_Merge_Key",
  all = FALSE,
  sort = FALSE,
  suffixes = c(
    "_Metadata",
    "_AST"
  )
)


# -----------------------------
# 11. Reorder important columns
# -----------------------------

priority_columns <- c(
  "Run_Accession",
  "BioProject",
  "BioSample_Merge_Key",
  "BioSample_Metadata",
  "BioSample_AST",
  "Sample_Name",
  "Country",
  "Collection_Year",
  "Collection_Date",
  "Isolation_Source",
  "Host_Category",
  "Host_Raw",
  "Host_Disease",
  "Organism",
  "Experiment",
  "LibraryLayout",
  "Platform",
  "Model",
  "Cohort_Status"
)

priority_columns <- priority_columns[
  priority_columns %in% names(master)
]

remaining_columns <- setdiff(
  names(master),
  priority_columns
)

master <- master[
  c(
    priority_columns,
    remaining_columns
  )
]


# -----------------------------
# 12. Sort master dataset
# -----------------------------

sort_columns <- intersect(
  c(
    "BioSample_Merge_Key",
    "Run_Accession",
    "Antibiotic",
    "Antimicrobial",
    "Drug"
  ),
  names(master)
)

if (length(sort_columns) > 0) {

  sorting_data <- master[
    sort_columns
  ]

  sorting_order <- do.call(
    order,
    c(
      sorting_data,
      list(
        na.last = TRUE
      )
    )
  )

  master <- master[
    sorting_order,
    ,
    drop = FALSE
  ]
}


# -----------------------------
# 13. Save master database
# -----------------------------

write.csv(
  master,
  master_output,
  row.names = FALSE,
  na = ""
)


# -----------------------------
# 14. Create unique run list
# -----------------------------

if (!is.null(run_column)) {

  run_list_columns <- intersect(
    c(
      run_column,
      metadata_biosample_column,
      "BioProject",
      "Sample_Name",
      "Country",
      "Collection_Year",
      "Isolation_Source",
      "Host_Category",
      "Organism",
      "LibraryLayout",
      "Platform",
      "Model",
      "Cohort_Status"
    ),
    names(metadata)
  )

  run_list <- unique(
    metadata[
      metadata$BioSample_Merge_Key %in%
        matched_biosamples,
      run_list_columns,
      drop = FALSE
    ]
  )

  run_list <- run_list[
    !is.na(run_list[[run_column]]) &
      trimws(run_list[[run_column]]) != "",
    ,
    drop = FALSE
  ]

} else {

  run_list <- unique(
    metadata[
      metadata$BioSample_Merge_Key %in%
        matched_biosamples,
      ,
      drop = FALSE
    ]
  )
}


write.csv(
  run_list,
  run_list_output,
  row.names = FALSE,
  na = ""
)


# -----------------------------
# 15. Calculate summary metrics
# -----------------------------

unique_master_biosamples <- length(
  unique(master$BioSample_Merge_Key)
)

unique_master_runs <- NA_integer_

if ("Run_Accession" %in% names(master)) {

  unique_master_runs <- length(
    unique(
      master$Run_Accession[
        !is.na(master$Run_Accession) &
          master$Run_Accession != ""
      ]
    )
  )
}


summary_table <- data.frame(

  Metric = c(
    "Human clinical metadata rows imported",
    "AST MIC rows imported",
    "Metadata rows missing BioSample",
    "AST rows missing BioSample",
    "Metadata exact duplicates removed",
    "AST exact duplicates removed",
    "Unique human metadata BioSamples",
    "Unique AST BioSamples",
    "Matched BioSamples",
    "Human metadata BioSamples without AST match",
    "AST BioSamples without human metadata match",
    "Master long-format rows",
    "Unique master BioSamples",
    "Unique master sequencing runs",
    "Run-list rows"
  ),

  Value = c(
    metadata_before_dedup,
    ast_before_dedup,
    metadata_missing_biosample,
    ast_missing_biosample,
    metadata_duplicates_removed,
    ast_duplicates_removed,
    length(metadata_biosamples),
    length(ast_biosamples),
    length(matched_biosamples),
    length(metadata_only_biosamples),
    length(ast_only_biosamples),
    nrow(master),
    unique_master_biosamples,
    unique_master_runs,
    nrow(run_list)
  ),

  stringsAsFactors = FALSE
)


write.csv(
  summary_table,
  summary_output,
  row.names = FALSE
)


# -----------------------------
# 16. Final validation
# -----------------------------

if (nrow(master) == 0) {

  stop(
    paste(
      "The merge produced zero rows.",
      "Check BioSample accession formatting in both inputs."
    )
  )
}

if (anyDuplicated(names(master)) > 0) {

  warning(
    "The master dataset contains duplicated column names."
  )
}


# -----------------------------
# 17. Print completion summary
# -----------------------------

cat("\n")
cat("=================================================\n")
cat("Module 06 completed successfully\n")
cat("=================================================\n")

cat(
  "Matched BioSamples:",
  length(matched_biosamples),
  "\n"
)

cat(
  "Master long-format rows:",
  nrow(master),
  "\n"
)

cat(
  "Unique master BioSamples:",
  unique_master_biosamples,
  "\n"
)

cat(
  "Unique sequencing runs:",
  unique_master_runs,
  "\n"
)

cat(
  "Human metadata BioSamples without AST match:",
  length(metadata_only_biosamples),
  "\n"
)

cat(
  "AST BioSamples without human metadata match:",
  length(ast_only_biosamples),
  "\n"
)

cat("\n")
cat("Outputs created:\n")
cat(" -", master_output, "\n")
cat(" -", run_list_output, "\n")
cat(" -", unmatched_metadata_output, "\n")
cat(" -", unmatched_ast_output, "\n")
cat(" -", summary_output, "\n")

