.libPaths(c("~/R/library", .libPaths()))

library(ape)
library(readr)
library(dplyr)

# Resolve repository root from script location
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
stopifnot(length(script_arg) == 1)
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))

TREE_FILE <- file.path(base_dir, "data", "KP_Final_1143_core_genome_tree.nwk")
MASTER_FILE <- file.path(base_dir, "data", "01_Stage4_Integrated_Isolate_Master_FINAL_1143.tsv")
GP_FILE <- file.path(base_dir, "data", "03_Isolate_Level_Siderophore_AMR_Genotype_Phenotype.tsv")
OUT <- file.path(base_dir, "results", "phylogeny", "stage2_amr_genotype_phenotype")

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

gp <- read_tsv(
  GP_FILE,
  show_col_types = FALSE
)

# ============================================================
# VALIDATION
# ============================================================

stopifnot(length(tree$tip.label) == 1143)
stopifnot(nrow(master) == 1143)
stopifnot(n_distinct(master$SRR) == 1143)

stopifnot(
  length(
    setdiff(master$SRR, tree$tip.label)
  ) == 0
)

cat("\n========================================\n")
cat("STAGE 2 PHYLOGENETIC ANALYSIS\n")
cat("SIDEROPHORE-AMR GENOTYPE-PHENOTYPE\n")
cat("========================================\n")

cat("Tree isolates:", length(tree$tip.label), "\n")
cat("Genomic isolates:", nrow(master), "\n")
cat(
  "Phenotype isolates:",
  n_distinct(gp$SRR),
  "\n"
)
cat(
  "Phenotype observations:",
  nrow(gp),
  "\n"
)

# ============================================================
# FULL TREE DISTANCE MATRIX
# ============================================================

D <- cophenetic.phylo(tree)

# ============================================================
# FUNCTION:
# PHYLOGENETIC CLUSTERING OF A POSITIVE GROUP
#
# Observed mean pairwise distance among positives is compared
# with random groups of identical size drawn from the tree.
#
# clustering_ratio < 1:
# positives are closer than expected.
# ============================================================

set.seed(20260812)

N_PERM <- 10000

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
        test_status = "Descriptive only: positive n below threshold"
      )
    )
  }

  obsD <- D[
    positive_ids,
    positive_ids,
    drop = FALSE
  ]

  obs_vals <- obsD[
    upper.tri(obsD)
  ]

  observed <- mean(
    obs_vals
  )

  perm_means <- numeric(
    N_PERM
  )

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

  expected <- mean(
    perm_means
  )

  # one-sided:
  # probability random sets are as close or closer
  p <- (
    sum(
      perm_means <= observed
    ) + 1
  ) / (
    N_PERM + 1
  )

  tibble(
    positive_n = n_pos,
    universe_n = length(universe_ids),
    observed_mean_distance = observed,
    expected_mean_distance = expected,
    clustering_ratio =
      observed / expected,
    permutation_P = p,
    test_status = "Tested"
  )
}

# ============================================================
# PART A
# AMR GENOTYPE PHYLOGENETIC STRUCTURE
# ============================================================

AMR_MARKERS <- c(
  "KPC",
  "NDM",
  "OXA48_like",
  "VIM_IMP",
  "CTXM",
  "qnr",
  "AAC6Ibcr",
  "gyrA_mut",
  "parC_mut",
  "mcr",
  "colistin_regulatory"
)

universe1143 <- tree$tip.label

amr_results <- list()

k <- 1

for(marker in AMR_MARKERS) {

  positive_ids <- master %>%
    filter(
      .data[[marker]] == TRUE
    ) %>%
    pull(SRR)

  z <- test_positive_clustering(
    positive_ids = positive_ids,
    universe_ids = universe1143,
    min_positive = 20
  )

  z$marker <- marker

  amr_results[[k]] <- z

  k <- k + 1
}

amr_phylo <- bind_rows(
  amr_results
) %>%
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
  amr_phylo,
  file.path(
    OUT,
    "01_AMR_Marker_Phylogenetic_Clustering.tsv"
  )
)

# ============================================================
# PART B
# GENOTYPE-PHENOTYPE DISCORDANCE PHYLOGENETIC STRUCTURE
#
# One record per isolate-antimicrobial already present.
# Test:
#   overall discordance
#   G+/P-S
#   G-/P-R
# ============================================================

gp2 <- gp %>%
  filter(
    Phenotype %in% c("R", "S")
  ) %>%
  distinct(
    SRR,
    Antibiotic,
    .keep_all = TRUE
  )

antibiotics <- c(
  "Ciprofloxacin",
  "Colistin",
  "Imipenem",
  "Meropenem"
)

discordance_results <- list()

k <- 1

for(drug in antibiotics) {

  z <- gp2 %>%
    filter(
      Antibiotic == drug
    )

  universe <- intersect(
    z$SRR,
    tree$tip.label
  )

  # ----------------------------------------------------------
  # Overall discordance
  # ----------------------------------------------------------

  ids_discordant <- z %>%
    filter(
      concordance == "Discordant"
    ) %>%
    pull(SRR)

  a <- test_positive_clustering(
    positive_ids = ids_discordant,
    universe_ids = universe,
    min_positive = 20
  )

  a$Antibiotic <- drug
  a$state <- "Overall discordance"

  discordance_results[[k]] <- a
  k <- k + 1

  # ----------------------------------------------------------
  # G+ / P-S
  # ----------------------------------------------------------

  ids_fp <- z %>%
    filter(
      integrated_state == "G+ / P-S"
    ) %>%
    pull(SRR)

  a <- test_positive_clustering(
    positive_ids = ids_fp,
    universe_ids = universe,
    min_positive = 20
  )

  a$Antibiotic <- drug
  a$state <- "G+ / P-S"

  discordance_results[[k]] <- a
  k <- k + 1

  # ----------------------------------------------------------
  # G- / P-R
  # ----------------------------------------------------------

  ids_fn <- z %>%
    filter(
      integrated_state == "G- / P-R"
    ) %>%
    pull(SRR)

  a <- test_positive_clustering(
    positive_ids = ids_fn,
    universe_ids = universe,
    min_positive = 20
  )

  a$Antibiotic <- drug
  a$state <- "G- / P-R"

  discordance_results[[k]] <- a
  k <- k + 1
}

discordance_phylo <- bind_rows(
  discordance_results
) %>%
  mutate(
    FDR_P = p.adjust(
      permutation_P,
      method = "BH"
    )
  ) %>%
  select(
    Antibiotic,
    state,
    everything()
  ) %>%
  arrange(
    is.na(FDR_P),
    FDR_P
  )

write_tsv(
  discordance_phylo,
  file.path(
    OUT,
    "02_Genotype_Phenotype_Discordance_Phylogenetic_Clustering.tsv"
  )
)

# ============================================================
# PART C
# SIDEROPHORE PROFILE DESCRIPTIVE PHYLOGENETIC CONTEXT
#
# This is descriptive/global context only.
# Formal Stage 1 within-ST testing has already been performed.
# ============================================================

siderophore_counts <- master %>%
  count(
    siderophore_profile,
    name = "isolates"
  ) %>%
  arrange(
    desc(isolates)
  )

write_tsv(
  siderophore_counts,
  file.path(
    OUT,
    "03_Siderophore_Profile_Counts_Global_Tree.tsv"
  )
)

# ============================================================
# SIGNIFICANT RESULTS
# ============================================================

sig_amr <- amr_phylo %>%
  filter(
    test_status == "Tested",
    !is.na(FDR_P),
    FDR_P < 0.05
  )

sig_discordance <- discordance_phylo %>%
  filter(
    test_status == "Tested",
    !is.na(FDR_P),
    FDR_P < 0.05
  )

write_tsv(
  sig_amr,
  file.path(
    OUT,
    "04_SIGNIFICANT_AMR_Phylogenetic_Clustering.tsv"
  )
)

write_tsv(
  sig_discordance,
  file.path(
    OUT,
    "05_SIGNIFICANT_Discordance_Phylogenetic_Clustering.tsv"
  )
)

# ============================================================
# REPORT
# ============================================================

cat("\n========================================\n")
cat("AMR MARKER PHYLOGENETIC STRUCTURE\n")
cat("========================================\n")

print(
  amr_phylo,
  n = Inf,
  width = Inf
)

cat("\n========================================\n")
cat("SIGNIFICANT AMR MARKERS\n")
cat("========================================\n")

if(nrow(sig_amr) == 0) {

  cat(
    "No AMR marker remained significant after BH-FDR correction.\n"
  )

} else {

  print(
    sig_amr,
    n = Inf,
    width = Inf
  )
}

cat("\n========================================\n")
cat("GENOTYPE-PHENOTYPE DISCORDANCE\n")
cat("========================================\n")

print(
  discordance_phylo,
  n = Inf,
  width = Inf
)

cat("\n========================================\n")
cat("SIGNIFICANT DISCORDANCE PATTERNS\n")
cat("========================================\n")

if(nrow(sig_discordance) == 0) {

  cat(
    "No discordance state remained significant after BH-FDR correction.\n"
  )

} else {

  print(
    sig_discordance,
    n = Inf,
    width = Inf
  )
}

cat("\n========================================\n")
cat("STAGE 2 PHYLOGENETIC ANALYSIS COMPLETE\n")
cat("========================================\n")


