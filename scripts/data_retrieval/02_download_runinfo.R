###############################################################
# Module 2: Download SRA RunInfo for all BioProjects
# Klebsiella pneumoniae Data Retrieval Pipeline
###############################################################

library(readr)
library(dplyr)
library(tibble)

#--------------------------------------------------------------
# File locations
#--------------------------------------------------------------

input_file <- "../outputs/01_bioproject_summary.csv"

output_file <- "../outputs/02_runinfo_master.csv"

log_file <- "../logs/02_runinfo_retrieval_log.csv"

checkpoint_file <- "../outputs/02_runinfo_checkpoint.csv"

#--------------------------------------------------------------
# Confirm Entrez Direct is available
#--------------------------------------------------------------

if (Sys.which("esearch") == "") {
  stop(
    paste(
      "Entrez Direct is not loaded.",
      "Run: module load StdEnv/2023 edirect"
    )
  )
}

if (Sys.which("efetch") == "") {
  stop(
    paste(
      "efetch is not available.",
      "Run: module load StdEnv/2023 edirect"
    )
  )
}

#--------------------------------------------------------------
# Read BioProject list
#--------------------------------------------------------------

bioproject_table <- read_csv(
  input_file,
  show_col_types = FALSE
)

if (!"BioProject" %in% names(bioproject_table)) {
  stop("The input file does not contain a BioProject column.")
}

bio_projects <- bioproject_table %>%
  filter(
    !is.na(BioProject),
    BioProject != ""
  ) %>%
  distinct(BioProject) %>%
  pull(BioProject)

cat("BioProjects found:", length(bio_projects), "\n")

#--------------------------------------------------------------
# Download RunInfo for one BioProject
#--------------------------------------------------------------

download_runinfo <- function(bioproject) {

  cat("Processing:", bioproject, "\n")

  command <- paste(
    "esearch -db sra -query",
    shQuote(bioproject),
    "| efetch -format runinfo"
  )

  runinfo_text <- tryCatch(
    {
      system(
        command,
        intern = TRUE,
        ignore.stderr = FALSE
      )
    },
    error = function(e) {
      structure(
        character(),
        error_message = conditionMessage(e)
      )
    }
  )

  error_message <- attr(runinfo_text, "error_message")

  if (
    length(runinfo_text) == 0 ||
    !startsWith(runinfo_text[1], "Run,")
  ) {

    return(
      list(
        data = NULL,
        log = tibble(
          BioProject = bioproject,
          Status = "No_RunInfo",
          Run_Count = 0L,
          Error = ifelse(
            is.null(error_message),
            NA_character_,
            error_message
          ),
          Retrieval_Time = as.character(Sys.time())
        )
      )
    )
  }

  runinfo_csv <- paste(
    runinfo_text,
    collapse = "\n"
  )

  runinfo <- tryCatch(
     {
    read_csv(
      I(runinfo_csv),
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "unique",
      col_types = cols(.default = col_character())
    )
    },
    error = function(e) {
      structure(
        NULL,
        error_message = conditionMessage(e)
      )
    }
  )

  if (is.null(runinfo)) {

    return(
      list(
        data = NULL,
        log = tibble(
          BioProject = bioproject,
          Status = "CSV_Parse_Failed",
          Run_Count = 0L,
          Error = attr(runinfo, "error_message"),
          Retrieval_Time = as.character(Sys.time())
        )
      )
    )
  }

  if (nrow(runinfo) == 0) {

    return(
      list(
        data = NULL,
        log = tibble(
          BioProject = bioproject,
          Status = "No_Runs",
          Run_Count = 0L,
          Error = NA_character_,
          Retrieval_Time = as.character(Sys.time())
        )
      )
    )
  }

  runinfo <- runinfo %>%
  mutate(
    across(everything(), as.character),
    Source_BioProject = as.character(bioproject),
    .before = 1
  )

  list(
    data = runinfo,
    log = tibble(
      BioProject = bioproject,
      Status = "Success",
      Run_Count = nrow(runinfo),
      Error = NA_character_,
      Retrieval_Time = as.character(Sys.time())
    )
  )
}

#--------------------------------------------------------------
# Process all BioProjects
#--------------------------------------------------------------

runinfo_results <- list()
retrieval_logs <- list()

for (i in seq_along(bio_projects)) {

  current_bioproject <- bio_projects[i]

  cat(
    "\n[",
    i,
    "/",
    length(bio_projects),
    "] ",
    current_bioproject,
    "\n",
    sep = ""
  )

  result <- download_runinfo(current_bioproject)

  runinfo_results[[i]] <- result$data
  retrieval_logs[[i]] <- result$log

  # Save progress after each BioProject
  current_logs <- bind_rows(retrieval_logs)

  write_csv(
    current_logs,
    log_file,
    na = ""
  )

  current_successful <- runinfo_results[
    !vapply(runinfo_results, is.null, logical(1))
  ]

  if (length(current_successful) > 0) {

    checkpoint_data <- bind_rows(current_successful)

    write_csv(
      checkpoint_data,
      checkpoint_file,
      na = ""
    )
  }

  Sys.sleep(0.4)
}

#--------------------------------------------------------------
# Merge successful RunInfo tables
#--------------------------------------------------------------

successful_results <- runinfo_results[
  !vapply(runinfo_results, is.null, logical(1))
]

retrieval_log <- bind_rows(retrieval_logs)

if (length(successful_results) == 0) {

  write_csv(
    retrieval_log,
    log_file,
    na = ""
  )

  stop("No RunInfo records were downloaded.")
}

runinfo_master <- bind_rows(successful_results)

if ("Run" %in% names(runinfo_master)) {

  runinfo_master <- runinfo_master %>%
    distinct(Run, .keep_all = TRUE)
}

#--------------------------------------------------------------
# Save final outputs
#--------------------------------------------------------------

write_csv(
  runinfo_master,
  output_file,
  na = ""
)

write_csv(
  retrieval_log,
  log_file,
  na = ""
)

# Remove temporary checkpoint after successful completion
if (file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
}

#--------------------------------------------------------------
# Final summary
#--------------------------------------------------------------

cat("\n========================================\n")
cat("Module 2 completed successfully\n")
cat("========================================\n")

cat(
  "BioProjects processed:",
  length(bio_projects),
  "\n"
)

cat(
  "Successful BioProjects:",
  sum(retrieval_log$Status == "Success"),
  "\n"
)

cat(
  "BioProjects without RunInfo:",
  sum(retrieval_log$Status != "Success"),
  "\n"
)

cat(
  "Unique SRA runs:",
  nrow(runinfo_master),
  "\n"
)

cat(
  "Master output:",
  output_file,
  "\n"
)

cat(
  "Retrieval log:",
  log_file,
  "\n"
)


