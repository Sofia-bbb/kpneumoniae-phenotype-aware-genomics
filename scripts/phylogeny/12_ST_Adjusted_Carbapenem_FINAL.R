

library(ape)
library(readr)
library(dplyr)
library(tibble)

set.seed(20260812)
N_PERM <- 10000

# Repository root
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))

TREE_FILE <- file.path(
  base_dir, "data", "KP_Final_1143_core_genome_tree.nwk"
)

GP_FILE <- file.path(
  base_dir, "data",
  "03_Isolate_Level_Siderophore_AMR_Genotype_Phenotype.tsv"
)

output_dir <- file.path(base_dir, "results", "phylogeny")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

OUT_FILE <- file.path(
  output_dir,
  "06_ST_Adjusted_Carbapenem_Gneg_PR_Clustering.tsv"
)

tree <- read.tree(TREE_FILE)
D <- cophenetic.phylo(tree)

gp <- read_tsv(
  GP_FILE,
  show_col_types = FALSE
) %>%
  filter(Phenotype %in% c("R","S")) %>%
  distinct(SRR, Antibiotic, .keep_all = TRUE) %>%
  filter(SRR %in% tree$tip.label)

test_ST_stratified <- function(drug) {

  z <- gp %>%
    filter(Antibiotic == drug)

  positives <- z %>%
    filter(integrated_state == "G- / P-R") %>%
    distinct(SRR, .keep_all = TRUE)

  positive_ids <- positives$SRR

  obsD <- D[
    positive_ids,
    positive_ids,
    drop = FALSE
  ]

  observed <- mean(
    obsD[upper.tri(obsD)]
  )

  st_counts <- positives %>%
    count(ST, name = "n_positive")

  perm_means <- numeric(N_PERM)

  for(i in seq_len(N_PERM)) {

    sampled_ids <- character(0)

    for(j in seq_len(nrow(st_counts))) {

      this_ST <- st_counts$ST[j]
      n_needed <- st_counts$n_positive[j]

      pool <- z %>%
        filter(ST == this_ST) %>%
        pull(SRR) %>%
        unique()

      if(length(pool) < n_needed) {
        stop(
          paste(
            "Insufficient isolates in",
            this_ST,
            "for",
            drug
          )
        )
      }

      sampled_ids <- c(
        sampled_ids,
        sample(
          pool,
          n_needed,
          replace = FALSE
        )
      )
    }

    PD <- D[
      sampled_ids,
      sampled_ids,
      drop = FALSE
    ]

    perm_means[i] <- mean(
      PD[upper.tri(PD)]
    )
  }

  expected <- mean(perm_means)

  p <- (
    sum(perm_means <= observed) + 1
  ) / (
    N_PERM + 1
  )

  tibble(
    Antibiotic = drug,
    positive_n = length(positive_ids),
    ST_count = n_distinct(positives$ST),
    observed_mean_distance = observed,
    ST_adjusted_expected_distance = expected,
    ST_adjusted_clustering_ratio = observed / expected,
    ST_adjusted_permutation_P = p
  )
}

results <- bind_rows(
  test_ST_stratified("Imipenem"),
  test_ST_stratified("Meropenem")
) %>%
  mutate(
    ST_adjusted_FDR_P = p.adjust(
      ST_adjusted_permutation_P,
      method = "BH"
    )
  )

write_tsv(
  results,
  OUT_FILE
)

cat("\n============================================\n")
cat("FINAL ST-STRATIFIED G-/P-R CLUSTERING\n")
cat("10,000 permutations\n")
cat("============================================\n\n")

print(results, n = Inf, width = Inf)

cat("\nSaved to:\n", OUT_FILE, "\n")



