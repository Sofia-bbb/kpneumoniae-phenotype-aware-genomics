#!/usr/bin/env Rscript

# ============================================================
# Module 4A: NCBI Pathogen Detection AST/MIC Retrieval
#
# Purpose:
#   Retrieve and standardize Klebsiella pneumoniae AST/MIC
#   records independently of Modules 1–3.
#
# Input:
#   inputs/ncbi_pathogen_ast.csv
#
# Outputs:
#   outputs/04_AST_MIC_long.csv
#   outputs/04_ast_biosamples.txt
#   logs/04_AST_MIC_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

ast_input <- "inputs/ncbi_pathogen_ast.csv"

ast_output <- "outputs/04_AST_MIC_long.csv"
biosample_output <- "outputs/04_ast_biosamples.txt"
summary_output <- "logs/04_AST_MIC_summary.csv"

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

if (!file.exists(ast_input)) {
  stop("Missing input file: ", ast_input)
}

message("Reading NCBI Pathogen Detection AST data...")

ast_raw <- read_csv(
  ast_input,
  show_col_types = FALSE,
  guess_max = 100000,
  progress = FALSE
)

message("Total AST records imported: ", nrow(ast_raw))

ast_clean <- ast_raw |>
  transmute(
    BioSample = str_to_upper(str_trim(`#BioSample`)),
    Scientific_Name = str_squish(`Scientific name`),
    Isolation_Type = str_squish(`Isolation type`),
    Location = str_squish(Location),
    Isolation_Source = str_squish(`Isolation source`),
    Isolate = str_squish(Isolate),
    Antibiotic_Original = str_squish(Antibiotic),
    Resistance_Phenotype_Original =
      str_squish(`Resistance phenotype`),
    Measurement_Sign = str_squish(`Measurement sign`),
    MIC_Value = suppressWarnings(as.numeric(`MIC (mg/L)`)),
    Disk_Diffusion_mm =
      suppressWarnings(as.numeric(`Disk diffusion (mm)`)),
    Testing_Platform =
      str_squish(`Laboratory typing platform`),
    Vendor = str_squish(Vendor),
    Method_or_Reagent =
      str_squish(
        `Laboratory typing method version or reagent`
      ),
    Testing_Standard = str_squish(`Testing standard`),
    Create_Date = str_squish(`Create date`)
  ) |>
  mutate(
    across(
      where(is.character),
      ~ na_if(.x, "")
    ),

    Antibiotic = case_when(
      str_to_lower(Antibiotic_Original) == "meropenem" ~
        "Meropenem",
      str_to_lower(Antibiotic_Original) == "imipenem" ~
        "Imipenem",
      str_to_lower(Antibiotic_Original) == "ciprofloxacin" ~
        "Ciprofloxacin",
      str_to_lower(Antibiotic_Original) == "colistin" ~
        "Colistin",
      TRUE ~ NA_character_
    ),

    Resistance_Phenotype = case_when(
      str_to_lower(Resistance_Phenotype_Original) %in%
        c("resistant", "r") ~ "R",

      str_to_lower(Resistance_Phenotype_Original) %in%
        c("intermediate", "i") ~ "I",

      str_to_lower(Resistance_Phenotype_Original) %in%
        c("susceptible", "sensitive", "s") ~ "S",

      TRUE ~ Resistance_Phenotype_Original
    ),

    MIC_Units = if_else(
      !is.na(MIC_Value),
      "mg/L",
      NA_character_
    ),

    MIC_Reported = case_when(
      !is.na(MIC_Value) &
        !is.na(Measurement_Sign) ~
        paste0(Measurement_Sign, MIC_Value),

      !is.na(MIC_Value) ~
        as.character(MIC_Value),

      TRUE ~ NA_character_
    )
  )

klebsiella_ast <- ast_clean |>
  filter(
    str_detect(
      str_to_lower(Scientific_Name),
      "^klebsiella pneumoniae"
    )
  )

target_ast <- klebsiella_ast |>
  filter(!is.na(Antibiotic)) |>
  filter(
    !is.na(Resistance_Phenotype) |
      !is.na(MIC_Value)
  ) |>
  distinct() |>
  arrange(
    BioSample,
    Antibiotic,
    MIC_Value
  )

write_csv(
  target_ast,
  ast_output
)

biosamples <- target_ast |>
  filter(!is.na(BioSample)) |>
  distinct(BioSample) |>
  arrange(BioSample)

write_lines(
  biosamples$BioSample,
  biosample_output
)

summary_table <- tibble(
  Metric = c(
    "Total AST records imported",
    "Klebsiella pneumoniae AST records",
    "Target AST or MIC records",
    "Unique BioSamples requiring SRA mapping",
    "Meropenem records",
    "Imipenem records",
    "Ciprofloxacin records",
    "Colistin records",
    "Records with MIC",
    "Records with phenotype",
    "Records with both MIC and phenotype"
  ),

  Value = c(
    nrow(ast_raw),
    nrow(klebsiella_ast),
    nrow(target_ast),
    n_distinct(target_ast$BioSample),

    sum(target_ast$Antibiotic == "Meropenem"),
    sum(target_ast$Antibiotic == "Imipenem"),
    sum(target_ast$Antibiotic == "Ciprofloxacin"),
    sum(target_ast$Antibiotic == "Colistin"),

    sum(!is.na(target_ast$MIC_Value)),
    sum(!is.na(target_ast$Resistance_Phenotype)),

    sum(
      !is.na(target_ast$MIC_Value) &
        !is.na(target_ast$Resistance_Phenotype)
    )
  )
)

write_csv(
  summary_table,
  summary_output
)

message("")
message("==============================================")
message("Module 4A completed successfully")
message("==============================================")
message("Target AST/MIC records: ", nrow(target_ast))
message("BioSamples for SRA mapping: ", nrow(biosamples))
message("")
message("Outputs:")
message("  ", ast_output)
message("  ", biosample_output)
message("")
message("Summary:")
message("  ", summary_output)
message("==============================================")

