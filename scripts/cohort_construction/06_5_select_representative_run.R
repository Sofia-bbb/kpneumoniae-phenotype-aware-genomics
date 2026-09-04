library(dplyr)
library(readr)
library(stringr)

# Read the human clinical sequencing-run list
runs <- read_csv(
  "outputs/06_Human_Clinical_Run_List.csv",
  show_col_types = FALSE
)

# Select one representative sequencing run per BioSample
selected <- runs %>%
  mutate(
    layout_priority = case_when(
      str_to_upper(LibraryLayout) == "PAIRED" ~ 1,
      TRUE ~ 2
    ),
    platform_priority = case_when(
      str_to_upper(Platform) == "ILLUMINA" ~ 1,
      TRUE ~ 2
    )
  ) %>%
  arrange(
    BioSample,
    layout_priority,
    platform_priority,
    Run_Accession
  ) %>%
  group_by(BioSample) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    -layout_priority,
    -platform_priority
  )

# Create output folders if needed
dir.create("outputs", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# Save the final selected genome-run list
write_csv(
  selected,
  "outputs/06_Selected_Genome_Run_List.csv"
)

# Create summary report
summary_table <- tibble(
  Metric = c(
    "Total candidate sequencing runs",
    "Unique BioSamples",
    "Selected representative runs",
    "Selected paired-end runs",
    "Selected Illumina runs"
  ),
  Value = c(
    nrow(runs),
    n_distinct(runs$BioSample),
    nrow(selected),
    sum(str_to_upper(selected$LibraryLayout) == "PAIRED", na.rm = TRUE),
    sum(str_to_upper(selected$Platform) == "ILLUMINA", na.rm = TRUE)
  )
)

write_csv(
  summary_table,
  "logs/06_5_selected_runs_summary.csv"
)

print(summary_table)

