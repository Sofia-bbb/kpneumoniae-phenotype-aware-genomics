#!/usr/bin/env Rscript

# ============================================================
# Module 3: BioSample Metadata Retrieval and Harmonization
#
# Pipeline:
# NCBI Pathogen Detection Genotype–Phenotype Pipeline
#
# Inputs:
#   outputs/02_SRA_RunInfo.csv
#   outputs/03_biosample_metadata.xml
#
# Outputs:
#   outputs/03_Metadata.csv
#   outputs/03_Metadata_BioSample_level.csv
#   logs/03_metadata_summary.csv
#   logs/03_missing_metadata.csv
# ============================================================

suppressPackageStartupMessages({
  library(xml2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

runinfo_file <- "outputs/02_SRA_RunInfo.csv"
xml_file <- "outputs/03_biosample_metadata.xml"

biosample_output <- "outputs/03_Metadata_BioSample_level.csv"
run_output <- "outputs/03_Metadata.csv"

summary_output <- "logs/03_metadata_summary.csv"
missing_output <- "logs/03_missing_metadata.csv"

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 2. Check inputs
# ------------------------------------------------------------

if (!file.exists(runinfo_file)) {
  stop("Missing RunInfo file: ", runinfo_file)
}

if (!file.exists(xml_file)) {
  stop("Missing BioSample XML file: ", xml_file)
}

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

normalize_missing <- function(x) {

  x <- str_squish(as.character(x))

  missing_terms <- c(
    "",
    "na",
    "n/a",
    "null",
    "none",
    "missing",
    "not available",
    "not applicable",
    "not collected",
    "not provided",
    "unknown"
  )

  x[str_to_lower(x) %in% missing_terms] <- NA_character_

  x
}

first_nonmissing <- function(x) {

  x <- normalize_missing(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  x[[1]]
}

find_attribute <- function(attribute_names,
                           attribute_values,
                           candidate_names) {

  normalized_names <- attribute_names |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_+|_+$", "")

  candidate_names <- candidate_names |>
    str_to_lower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("^_+|_+$", "")

  matches <- which(normalized_names %in% candidate_names)

  if (length(matches) == 0) {
    return(NA_character_)
  }

  first_nonmissing(attribute_values[matches])
}

extract_year <- function(x) {

  x <- normalize_missing(x)

  year <- str_extract(
    x,
    "(?<![0-9])(19|20)[0-9]{2}(?![0-9])"
  )

  suppressWarnings(as.integer(year))
}

clean_country <- function(x) {

  x <- normalize_missing(x)

  x <- str_replace(
    x,
    ":.*$",
    ""
  )

  x <- str_squish(x)

  case_when(
    str_detect(str_to_lower(x), "^usa$|^united states") ~
      "United States",

    str_detect(str_to_lower(x), "^uk$|^united kingdom") ~
      "United Kingdom",

    str_detect(str_to_lower(x), "^canada$") ~
      "Canada",

    TRUE ~ x
  )
}

classify_host <- function(x) {

  x_lower <- str_to_lower(normalize_missing(x))

  case_when(
    is.na(x_lower) ~ NA_character_,

    str_detect(
      x_lower,
      "homo sapiens|human|patient|clinical"
    ) ~ "Human",

    TRUE ~ "Non-human_or_other"
  )
}

# ------------------------------------------------------------
# 4. Read SRA RunInfo
# ------------------------------------------------------------

message("Reading SRA RunInfo...")

runinfo <- read_csv(
  runinfo_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  mutate(
    BioSample = str_to_upper(
      normalize_missing(BioSample)
    ),
    BioProject = str_to_upper(
      normalize_missing(BioProject)
    ),
    Run = str_to_upper(
      normalize_missing(Run)
    )
  )

message("SRA runs imported: ", nrow(runinfo))
message(
  "Unique BioSamples in RunInfo: ",
  n_distinct(runinfo$BioSample)
)

# ------------------------------------------------------------
# 5. Read BioSample XML
# ------------------------------------------------------------

message("Reading BioSample XML...")

xml_document <- read_xml(
  xml_file,
  options = c("RECOVER", "HUGE")
)

biosample_nodes <- xml_find_all(
  xml_document,
  ".//BioSample"
)

message(
  "BioSample XML records found: ",
  length(biosample_nodes)
)

# ------------------------------------------------------------
# 6. Extract metadata from each BioSample
# ------------------------------------------------------------

message("Extracting BioSample metadata...")

biosample_metadata <- map_dfr(
  biosample_nodes,

  function(node) {

    accession <- xml_attr(node, "accession")
    sample_name <- xml_attr(node, "id")

    description_node <- xml_find_first(
      node,
      ".//Description/Title"
    )

    description_title <- if (
      inherits(description_node, "xml_missing")
    ) {
      NA_character_
    } else {
      xml_text(description_node)
    }

    organism_node <- xml_find_first(
      node,
      ".//Organism"
    )

    organism_name <- if (
      inherits(organism_node, "xml_missing")
    ) {
      NA_character_
    } else {
      xml_attr(organism_node, "taxonomy_name")
    }

    attribute_nodes <- xml_find_all(
      node,
      ".//Attributes/Attribute"
    )

    attribute_names <- map_chr(
      attribute_nodes,
      function(attribute_node) {

        name <- xml_attr(
          attribute_node,
          "attribute_name"
        )

        harmonized_name <- xml_attr(
          attribute_node,
          "harmonized_name"
        )

        if (
          is.na(name) ||
          name == ""
        ) {
          name <- harmonized_name
        }

        name
      }
    )

    attribute_values <- map_chr(
      attribute_nodes,
      xml_text
    )

    sample_name_attribute <- find_attribute(
      attribute_names,
      attribute_values,
      c(
        "sample_name",
        "sample name",
        "isolate",
        "strain",
        "isolate_name",
        "isolate name"
      )
    )

    country_raw <- find_attribute(
      attribute_names,
      attribute_values,
      c(
        "geo_loc_name",
        "geographic location",
        "geographic_location",
        "country",
        "location"
      )
    )

    collection_date_raw <- find_attribute(
      attribute_names,
      attribute_values,
      c(
        "collection_date",
        "collection date",
        "isolation_date",
        "isolation date",
        "collection_year",
        "year"
      )
    )

    isolation_source <- find_attribute(
      attribute_names,
      attribute_values,
      c(
        "isolation_source",
        "isolation source",
        "sample_source",
        "sample source",
        "specimen",
        "specimen_source",
        "specimen source",
        "source"
      )
    )

    host_raw <- find_attribute(
      attribute_names,
      attribute_values,
      c(
        "host",
        "host_subject_id",
        "host scientific name",
        "host_scientific_name",
        "host_common_name",
        "host common name"
      )
    )

    host_disease <- find_attribute(
      attribute_names,
      attribute_values,
      c(
        "host_disease",
        "host disease",
        "disease",
        "health_state",
        "health state"
      )
    )

    tibble(
      BioSample = str_to_upper(
        normalize_missing(accession)
      ),

      Sample_Name = first_nonmissing(
        c(
          sample_name_attribute,
          sample_name,
          description_title
        )
      ),

      Organism = normalize_missing(
        organism_name
      ),

      Country_Raw = normalize_missing(
        country_raw
      ),

      Country = clean_country(
        country_raw
      ),

      Collection_Date = normalize_missing(
        collection_date_raw
      ),

      Collection_Year = extract_year(
        collection_date_raw
      ),

      Isolation_Source = normalize_missing(
        isolation_source
      ),

      Host_Raw = normalize_missing(
        host_raw
      ),

      Host_Category = classify_host(
        host_raw
      ),

      Host_Disease = normalize_missing(
        host_disease
      )
    )
  }
)

biosample_metadata <- biosample_metadata |>
  distinct(BioSample, .keep_all = TRUE) |>
  arrange(BioSample)

write_csv(
  biosample_metadata,
  biosample_output
)

# ------------------------------------------------------------
# 7. Join metadata with RunInfo
# ------------------------------------------------------------

metadata_run_level <- runinfo |>
  left_join(
    biosample_metadata,
    by = "BioSample"
  ) |>
  transmute(
    Run_Accession = Run,
    BioProject,
    BioSample,
    Sample_Name,
    Country,
    Collection_Year,
    Isolation_Source,
    Host_Category,
    Host_Raw,
    Host_Disease,
    Organism,
    Collection_Date,
    Experiment,
    LibraryLayout,
    Platform,
    Model
  ) |>
  distinct() |>
  arrange(
    BioProject,
    BioSample,
    Run_Accession
  )

write_csv(
  metadata_run_level,
  run_output
)

# ------------------------------------------------------------
# 8. Missing metadata report
# ------------------------------------------------------------

missing_metadata <- metadata_run_level |>
  mutate(
    Missing_Country = is.na(Country),
    Missing_Collection_Year = is.na(Collection_Year),
    Missing_Isolation_Source = is.na(Isolation_Source),
    Missing_Host = is.na(Host_Category)
  ) |>
  filter(
    Missing_Country |
      Missing_Collection_Year |
      Missing_Isolation_Source |
      Missing_Host
  )

write_csv(
  missing_metadata,
  missing_output
)

# ------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------

summary_table <- tibble(
  Metric = c(
    "SRA runs imported",
    "Unique BioSamples in RunInfo",
    "BioSample XML records",
    "BioSamples parsed",
    "Run-level metadata records",
    "Unique runs in final metadata",
    "Unique BioSamples in final metadata",
    "Unique BioProjects",
    "Records with country",
    "Records with collection year",
    "Records with isolation source",
    "Records classified as human",
    "Records classified as non-human or other",
    "Records with missing host classification",
    "Canadian records"
  ),

  Value = c(
    nrow(runinfo),
    n_distinct(runinfo$BioSample),
    length(biosample_nodes),
    nrow(biosample_metadata),
    nrow(metadata_run_level),
    n_distinct(metadata_run_level$Run_Accession),
    n_distinct(metadata_run_level$BioSample),
    n_distinct(metadata_run_level$BioProject),
    sum(!is.na(metadata_run_level$Country)),
    sum(!is.na(metadata_run_level$Collection_Year)),
    sum(!is.na(metadata_run_level$Isolation_Source)),
    sum(metadata_run_level$Host_Category == "Human",
        na.rm = TRUE),
    sum(
      metadata_run_level$Host_Category ==
        "Non-human_or_other",
      na.rm = TRUE
    ),
    sum(is.na(metadata_run_level$Host_Category)),
    sum(metadata_run_level$Country == "Canada",
        na.rm = TRUE)
  )
)

write_csv(
  summary_table,
  summary_output
)

# ------------------------------------------------------------
# 10. Completion report
# ------------------------------------------------------------

message("")
message("==============================================")
message("Module 3 completed successfully")
message("==============================================")
message(
  "BioSample metadata records: ",
  nrow(biosample_metadata)
)
message(
  "Run-level metadata records: ",
  nrow(metadata_run_level)
)
message("")
message("Outputs:")
message("  ", biosample_output)
message("  ", run_output)
message("")
message("Logs:")
message("  ", summary_output)
message("  ", missing_output)
message("==============================================")

