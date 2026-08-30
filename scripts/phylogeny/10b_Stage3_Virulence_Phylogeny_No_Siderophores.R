library(ape)
library(readr)
library(dplyr)

# Repository root
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))

TREE_FILE <- file.path(
  base_dir, "data", "KP_Final_1143_core_genome_tree.nwk"
)

MASTER_FILE <- file.path(
  base_dir, "data",
  "01_Stage4_Integrated_Isolate_Master_FINAL_1143.tsv"
)

CORRECTED_VF_FILE <- file.path(
  base_dir, "data",
  "05_Isolate_Level_VF_Counts_No_Siderophores.tsv"
)

OUT <- file.path(
  base_dir, "results", "phylogeny",
  "No_Siderophore_VF_Burden"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

tree <- read.tree(TREE_FILE)

master <- read_tsv(
  MASTER_FILE,
  show_col_types = FALSE
)

corrected_vf <- read_tsv(
  CORRECTED_VF_FILE,
  show_col_types = FALSE
) %>%
  select(
    SRR,
    vf_factor_count_no_siderophores
  )

master <- master %>%
  left_join(
    corrected_vf,
    by = "SRR"
  )

# ============================================================
# VALIDATION
# ============================================================

stopifnot(length(tree$tip.label) == 1143)
stopifnot(nrow(master) == 1143)
stopifnot(n_distinct(master$SRR) == 1143)
stopifnot(
  sum(is.na(master$vf_factor_count_no_siderophores)) == 0
)

stopifnot(
  length(setdiff(master$SRR, tree$tip.label)) == 0
)

stopifnot(
  length(setdiff(tree$tip.label, master$SRR)) == 0
)

cat("\n========================================\n")
cat("VIRULENCE PHYLOGENETIC ANALYSIS\n")
cat("NON-SIDEROPHORE VF BURDEN\n")
cat("========================================\n")

cat("Tree isolates:", length(tree$tip.label), "\n")
cat("Master isolates:", nrow(master), "\n")

# ============================================================
# NON-SIDEROPHORE VIRULENCE BURDEN THRESHOLD
# ============================================================

vf_q3 <- quantile(
  master$vf_factor_count_no_siderophores,
  probs = 0.75,
  na.rm = TRUE,
  names = FALSE
)

master <- master %>%
  mutate(
    high_non_siderophore_VF_burden =
      vf_factor_count_no_siderophores >= vf_q3
  )

cat(
  "Non-siderophore VF burden Q3 threshold:",
  vf_q3,
  "\n"
)

cat(
  "High non-siderophore VF burden isolates:",
  sum(
    master$high_non_siderophore_VF_burden,
    na.rm = TRUE
  ),
  "\n"
)

# ============================================================
# FULL TREE DISTANCE MATRIX
# ============================================================

D <- cophenetic.phylo(tree)

set.seed(20260812)

N_PERM <- 10000

# ============================================================
# CLUSTERING FUNCTION
# ============================================================

test_positive_clustering <- function(
  positive_ids,
  universe_ids,
  min_positive = 20
) {

  positive_ids <- intersect(
    unique(positive_ids),
    universe_ids
  )

  n_pos <- length(positive_ids)

  if(n_pos < min_positive) {

    return(
      tibble(
        positive_n = n_pos,
        universe_n = length(universe_ids),
        observed_mean_distance = NA_real_,
        expected_mean_distance = NA_real_,
        clustering_ratio = NA_real_,
        permutation_P = NA_real_,
        test_status =
          "Descriptive only: positive n below threshold"
      )
    )
  }

  obsD <- D[
    positive_ids,
    positive_ids,
    drop = FALSE
  ]

  observed <- mean(
    obsD[
      upper.tri(obsD)
    ]
  )

  perm_means <- numeric(N_PERM)

  for(i in seq_len(N_PERM)) {

    ids <- sample(
      universe_ids,
      n_pos,
      replace = FALSE
    )

    PD <- D[
      ids,
      ids,
      drop = FALSE
    ]

    perm_means[i] <- mean(
      PD[
        upper.tri(PD)
      ]
    )
  }

  expected <- mean(perm_means)

  p <- (
    sum(perm_means <= observed) + 1
  ) / (
    N_PERM + 1
  )

  tibble(
    positive_n = n_pos,
    universe_n = length(universe_ids),
    observed_mean_distance = observed,
    expected_mean_distance = expected,
    clustering_ratio = observed / expected,
    permutation_P = p,
    test_status = "Tested"
  )
}

# ============================================================
# VIRULENCE FEATURES
# ============================================================

markers <- c(
  "colibactin_positive",
  "RmpADC_positive",
  "rmpA2_positive",
  "high_non_siderophore_VF_burden"
)

universe <- tree$tip.label

results <- list()
k <- 1

for(marker in markers) {

  positive_ids <- master %>%
    filter(.data[[marker]] == TRUE) %>%
    pull(SRR)

  z <- test_positive_clustering(
    positive_ids = positive_ids,
    universe_ids = universe,
    min_positive = 20
  )

  z$marker <- marker

  results[[k]] <- z
  k <- k + 1
}

virulence_phylo <- bind_rows(results) %>%
  mutate(
    FDR_P = p.adjust(
      permutation_P,
      method = "BH"
    )
  ) %>%
  select(
    marker,
    everything()
  ) %>%
  arrange(
    is.na(FDR_P),
    FDR_P
  )

write_tsv(
  virulence_phylo,
  file.path(
    OUT,
    "01_Virulence_Phylogenetic_Clustering_No_Siderophore_VF.tsv"
  )
)

# ============================================================
# PROFILE SUMMARY
# ============================================================

profile_summary <- master %>%
  group_by(siderophore_profile) %>%
  summarise(
    isolates = n(),

    colibactin_n =
      sum(colibactin_positive, na.rm = TRUE),

    RmpADC_n =
      sum(RmpADC_positive, na.rm = TRUE),

    rmpA2_n =
      sum(rmpA2_positive, na.rm = TRUE),

    high_non_siderophore_VF_n =
      sum(
        high_non_siderophore_VF_burden,
        na.rm = TRUE
      ),

    median_non_siderophore_VF =
      median(
        vf_factor_count_no_siderophores,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(desc(isolates))

write_tsv(
  profile_summary,
  file.path(
    OUT,
    "02_Virulence_By_Siderophore_Profile_No_Siderophore_VF.tsv"
  )
)

cat("\n========================================\n")
cat("PHYLOGENETIC RESULTS\n")
cat("========================================\n")

print(
  virulence_phylo,
  n = Inf,
  width = Inf
)

cat("\n========================================\n")
cat("PROFILE SUMMARY\n")
cat("========================================\n")

print(
  profile_summary,
  n = Inf,
  width = Inf
)


