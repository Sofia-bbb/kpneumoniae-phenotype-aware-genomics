#!/usr/bin/env Rscript

# ============================================================
# Module 05: Human Clinical Cohort Selection
# NCBI Pathogen Detection Genotype-Phenotype Pipeline
#
# Purpose:
#   Classify metadata records as:
#     1. Human_clinical
#     2. Excluded_nonhuman
#     3. Manual_review
#
# Input:
#   outputs/03_Metadata.csv
#
# Outputs:
#   outputs/05_Human_Clinical.csv
#   outputs/05_Excluded_NonHuman.csv
#   outputs/05_Manual_Review.csv
#   outputs/05_All_Classified.csv
#   logs/05_human_clinical_summary.csv
# ============================================================


# -----------------------------
# 1. File locations
# -----------------------------

input_file <- "outputs/03_Metadata.csv"

human_output <- "outputs/05_Human_Clinical.csv"
excluded_output <- "outputs/05_Excluded_NonHuman.csv"
review_output <- "outputs/05_Manual_Review.csv"
classified_output <- "outputs/05_All_Classified.csv"
summary_output <- "logs/05_human_clinical_summary.csv"


# -----------------------------
# 2. Check input file
# -----------------------------

if (!file.exists(input_file)) {
  stop(
    paste(
      "Input file was not found:",
      input_file
    )
  )
}


# -----------------------------
# 3. Read metadata
# -----------------------------

metadata <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A", "na", "n/a", "null", "NULL")
)

cat("Metadata rows imported:", nrow(metadata), "\n")
cat("Metadata columns imported:", ncol(metadata), "\n")


# -----------------------------
# 4. Normalize column names
# -----------------------------

normalize_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  return(x)
}

normalized_columns <- normalize_name(names(metadata))


# -----------------------------
# 5. Find relevant columns
# -----------------------------

find_column <- function(possible_names) {

  possible_names <- normalize_name(possible_names)

  match_index <- match(
    possible_names,
    normalized_columns,
    nomatch = 0
  )

  match_index <- match_index[match_index > 0]

  if (length(match_index) == 0) {
    return(NULL)
  }

  names(metadata)[match_index[1]]
}


host_column <- find_column(
  c(
    "Host_Category",
    "Host",
    "host_name",
    "host_scientific_name"
  )
)

host_raw_column <- find_column(
  c(
    "Host_Raw",
    "host_raw",
    "host_scientific_name"
  )
)

source_column <- find_column(
  c(
    "Isolation_Source",
    "isolation source",
    "source",
    "sample_source"
  )
)

sample_column <- find_column(
  c(
    "Sample_Name",
    "sample name",
    "sample_title",
    "title"
  )
)

host_disease_column <- find_column(
  c(
    "Host_Disease",
    "host disease",
    "disease"
  )
)


cat(
  "Host column:",
  ifelse(is.null(host_column), "Not found", host_column),
  "\n"
)

cat(
  "Isolation-source column:",
  ifelse(is.null(source_column), "Not found", source_column),
  "\n"
)

cat(
  "Sample-name column:",
  ifelse(is.null(sample_column), "Not found", sample_column),
  "\n"
)

cat(
  "Host-disease column:",
  ifelse(
    is.null(host_disease_column),
    "Not found",
    host_disease_column
  ),
  "\n"
)


# -----------------------------
# 6. Safely retrieve text fields
# -----------------------------

get_text_column <- function(column_name) {

  if (is.null(column_name)) {
    return(rep("", nrow(metadata)))
  }

  values <- metadata[[column_name]]
  values[is.na(values)] <- ""

  return(as.character(values))
}


host_text <- get_text_column(host_column)
host_raw_text <- get_text_column(host_raw_column)
source_text <- get_text_column(source_column)
sample_text <- get_text_column(sample_column)
disease_text <- get_text_column(host_disease_column)


# Combine fields for classification

classification_text <- paste(
  host_text,
  host_raw_text,
  source_text,
  sample_text,
  disease_text,
  sep = " | "
)

classification_text <- tolower(classification_text)


# -----------------------------
# 7. Human-clinical keywords
# -----------------------------

human_patterns <- paste(
  c(
    "homo sapiens",
    "\\bhuman\\b",
    "\\bpatient\\b",
    "clinical isolate",
    "\\bclinical\\b",
    "\\bblood\\b",
    "bloodstream",
    "\\bserum\\b",
    "\\bplasma\\b",
    "\\burine\\b",
    "urinary",
    "\\bsputum\\b",
    "respiratory",
    "bronchial",
    "bronchoalveolar",
    "\\bbal\\b",
    "tracheal",
    "nasopharyngeal",
    "oropharyngeal",
    "\\bthroat\\b",
    "\\bwound\\b",
    "\\babscess\\b",
    "\\bpus\\b",
    "cerebrospinal",
    "\\bcsf\\b",
    "\\bstool\\b",
    "\\bfeces\\b",
    "\\bfaeces\\b",
    "\\brectal\\b",
    "rectum",
    "perianal",
    "perirectal",
    "\\bgut\\b",
    "intestinal",
    "\\bcolon\\b",
    "peritoneal",
    "\\bascites\\b",
    "pleural",
    "synovial",
    "vaginal",
    "cervical",
    "urethral",
    "catheter",
    "surgical",
    "biopsy",
    "\\btissue\\b",
    "hospitalized",
    "hospitalised",
    "infection",
    "infected",
"hospital",
"hospitalized",
"hospitalised",
"clinic",
"clinical",
"healthcare",
"medical center",
"medical centre",
"inpatient",
"outpatient",
"icu",
"intensive care",
"ward",
"emergency department"
  ),
  collapse = "|"
)


# -----------------------------
# 8. Non-human exclusion keywords
# -----------------------------

nonhuman_patterns <- paste(
  c(
    "\\bfood\\b",
    "food product",
    "food sample",
    "ready.to.eat",
    "\\brte\\b",
    "retail food",
    "coconut",
    "\\bmeat\\b",
    "\\bbeef\\b",
    "\\bpork\\b",
    "\\bchicken\\b",
    "\\bturkey\\b",
    "\\bseafood\\b",
    "\\bfish\\b",
    "\\bshrimp\\b",
    "\\bmilk\\b",
    "\\bcheese\\b",
    "\\begg\\b",
    "\\bvegetable\\b",
    "\\bfruit\\b",
    "\\banimal\\b",
    "veterinary",
    "\\bbovine\\b",
    "\\bcattle\\b",
    "\\bovine\\b",
    "\\bswine\\b",
    "\\bpig\\b",
    "\\bporcine\\b",
    "\\bpoultry\\b",
    "\\bavian\\b",
    "\\bcanine\\b",
    "\\bdog\\b",
    "\\bfeline\\b",
    "\\bcat\\b",
    "\\bhorse\\b",
    "\\bequine\\b",
    "\\bgoat\\b",
    "\\bsheep\\b",
    "\\bmouse\\b",
    "\\bmice\\b",
    "\\brat\\b",
    "\\bbird\\b",
    "\\bwildlife\\b",
    "\\bwater\\b",
    "wastewater",
    "sewage",
    "effluent",
    "\\briver\\b",
    "\\blake\\b",
    "\\bsoil\\b",
    "\\bsediment\\b",
    "\\bplant\\b",
    "environmental",
    "environment",
    "retail environment",
    "surface swab",
    "facility surface",
    "equipment surface"
  ),
  collapse = "|"
)


# -----------------------------
# 9. Apply classification rules
# -----------------------------

human_match <- grepl(
  human_patterns,
  classification_text,
  perl = TRUE
)

nonhuman_match <- grepl(
  nonhuman_patterns,
  classification_text,
  perl = TRUE
)


metadata$Human_Evidence <- human_match
metadata$NonHuman_Evidence <- nonhuman_match

metadata$Cohort_Status <- "Manual_review"


# Explicit non-human evidence has priority

metadata$Cohort_Status[
  nonhuman_match
] <- "Excluded_nonhuman"


# Retain records with human evidence and no non-human evidence

metadata$Cohort_Status[
  human_match & !nonhuman_match
] <- "Human_clinical"


# Records containing both types of evidence require review

metadata$Cohort_Status[
  human_match & nonhuman_match
] <- "Manual_review"


# Store the exact text used for classification

metadata$Classification_Text <- classification_text


# -----------------------------
# 10. Split the datasets
# -----------------------------

human_clinical <- metadata[
  metadata$Cohort_Status == "Human_clinical",
  ,
  drop = FALSE
]

excluded_nonhuman <- metadata[
  metadata$Cohort_Status == "Excluded_nonhuman",
  ,
  drop = FALSE
]

manual_review <- metadata[
  metadata$Cohort_Status == "Manual_review",
  ,
  drop = FALSE
]


# -----------------------------
# 11. Save outputs
# -----------------------------

write.csv(
  human_clinical,
  human_output,
  row.names = FALSE,
  na = ""
)

write.csv(
  excluded_nonhuman,
  excluded_output,
  row.names = FALSE,
  na = ""
)

write.csv(
  manual_review,
  review_output,
  row.names = FALSE,
  na = ""
)

write.csv(
  metadata,
  classified_output,
  row.names = FALSE,
  na = ""
)


# -----------------------------
# 12. Create summary report
# -----------------------------

summary_table <- data.frame(
  Metric = c(
    "Total metadata rows",
    "Human clinical rows",
    "Excluded non-human rows",
    "Manual-review rows",
    "Human clinical percentage"
  ),

  Value = c(
    nrow(metadata),
    nrow(human_clinical),
    nrow(excluded_nonhuman),
    nrow(manual_review),
    round(
      100 * nrow(human_clinical) / nrow(metadata),
      2
    )
  ),

  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  summary_output,
  row.names = FALSE
)


# -----------------------------
# 13. Print completion summary
# -----------------------------

cat("\n")
cat("============================================\n")
cat("Module 05 completed successfully\n")
cat("============================================\n")
cat("Total rows:", nrow(metadata), "\n")
cat("Human clinical:", nrow(human_clinical), "\n")
cat("Excluded non-human:", nrow(excluded_nonhuman), "\n")
cat("Manual review:", nrow(manual_review), "\n")
cat("\n")
cat("Outputs created:\n")
cat(" -", human_output, "\n")
cat(" -", excluded_output, "\n")
cat(" -", review_output, "\n")
cat(" -", classified_output, "\n")
cat(" -", summary_output, "\n")




