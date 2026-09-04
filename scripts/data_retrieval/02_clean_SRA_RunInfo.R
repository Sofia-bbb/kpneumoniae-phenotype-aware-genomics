#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

input_file  <- "temp/04B_runinfo_raw.csv"
output_file <- "outputs/02_SRA_RunInfo.csv"
summary_file <- "logs/02_SRA_RunInfo_summary.csv"

dir.create("outputs", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

cat("Reading SRA RunInfo...\n")

runinfo <- read_csv(
  input_file,
  show_col_types = FALSE
)

cat("Records:", nrow(runinfo), "\n")

clean <- runinfo %>%
  transmute(
    Run = Run,
    BioProject = BioProject,
    BioSample = BioSample,
    Experiment = Experiment,
    LibraryLayout = LibraryLayout,
    Platform = Platform,
    Model = Model
  ) %>%
  distinct() %>%
  arrange(BioSample, Run)

write_csv(clean, output_file)

summary <- tibble(
  Metric = c(
    "Raw records",
    "Unique runs",
    "Unique BioSamples",
    "Unique BioProjects"
  ),
  Value = c(
    nrow(runinfo),
    n_distinct(clean$Run),
    n_distinct(clean$BioSample),
    n_distinct(clean$BioProject)
  )
)

write_csv(summary, summary_file)

cat("\n=========================================\n")
cat("SRA RunInfo cleaning completed\n")
cat("=========================================\n")
cat("Unique Runs:", n_distinct(clean$Run), "\n")
cat("Unique BioSamples:", n_distinct(clean$BioSample), "\n")
cat("Unique BioProjects:", n_distinct(clean$BioProject), "\n")
cat("\nOutput:\n")
cat(output_file, "\n")
cat(summary_file, "\n")

