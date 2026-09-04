#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)

work_dir <- Sys.getenv("KP_WORKDIR", unset = file.path(getwd(), "work"))
quast_dir <- file.path(work_dir, "quast")
output_file <- file.path(getwd(), "results", "Supplementary_Table_S2_Assembly_QC.csv")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

reports <- list.files(
  quast_dir,
  pattern = "report.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

extract_metric <- function(df, metric) {
  value <- df[df$Metric == metric, "Value", drop = TRUE]
  if (length(value) == 0) return(NA_character_)
  value[1]
}

read_quast <- function(file) {

  srr <- basename(dirname(file))

  x <- read.delim(
    file,
    header = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )

  x <- x[, 1:2]
  colnames(x) <- c("Metric", "Value")

  tibble(
    SRA_Accession = srr,
    Contigs = as.numeric(extract_metric(x, "# contigs")),
    Contigs_ge_1000bp = as.numeric(extract_metric(x, "# contigs (>= 1000 bp)")),
    Largest_Contig_bp = as.numeric(extract_metric(x, "Largest contig")),
    Total_Length_bp = as.numeric(extract_metric(x, "Total length")),
    GC_percent = as.numeric(extract_metric(x, "GC (%)")),
    N50_bp = as.numeric(extract_metric(x, "N50")),
    N75_bp = as.numeric(extract_metric(x, "N75")),
    L50 = as.numeric(extract_metric(x, "L50")),
    L75 = as.numeric(extract_metric(x, "L75")),
    Ns_per_100kbp = as.numeric(extract_metric(x, "# N's per 100 kbp"))
  )
}

summary <- map_dfr(reports, read_quast) %>%
  arrange(SRA_Accession)

write_csv(summary, output_file)

cat("QUAST reports processed:", nrow(summary), "\n")
cat("Output:", output_file, "\n")
