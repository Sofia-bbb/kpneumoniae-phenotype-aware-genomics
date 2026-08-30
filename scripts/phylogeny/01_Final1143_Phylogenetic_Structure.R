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

ELIGIBLE_FILE <- file.path(
  base_dir, "data",
  "02_Eligible_ST_Siderophore_Groups_nGE5.tsv"
)

OUT <- file.path(base_dir, "results", "phylogeny")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD + VALIDATE FINAL 1143 TREE
# ============================================================

tree <- read.tree(TREE_FILE)

master <- read_tsv(
  MASTER_FILE,
  show_col_types=FALSE
)

stopifnot(length(tree$tip.label) == 1143)
stopifnot(n_distinct(tree$tip.label) == 1143)
stopifnot(nrow(master) == 1143)
stopifnot(n_distinct(master$SRR) == 1143)

tree_ids <- sort(tree$tip.label)
master_ids <- sort(master$SRR)

missing_tree <- setdiff(master_ids, tree_ids)
extra_tree   <- setdiff(tree_ids, master_ids)

cat("\n========================================\n")
cat("FINAL 1143 PHYLOGENETIC INTEGRATION\n")
cat("========================================\n")
cat("Tree tips:", length(tree_ids), "\n")
cat("Master isolates:", length(master_ids), "\n")
cat("Master isolates missing from tree:", length(missing_tree), "\n")
cat("Tree tips absent from master:", length(extra_tree), "\n")

if(length(missing_tree) > 0 || length(extra_tree) > 0){
  stop("Tree/master mismatch detected.")
}

# ============================================================
# 2. WRITE FINAL TREE METADATA
# ============================================================

tree_meta <- master %>%
  select(
    SRR,
    ST,
    siderophore_profile,
    virulence_score,
    resistance_score,
    num_resistance_genes,
    amr_gene_count,
    KPC,
    NDM,
    OXA48_like,
    CTXM,
    qnr,
    AAC6Ibcr,
    gyrA_mut,
    parC_mut,
    mcr,
    colistin_regulatory,
    colibactin_positive,
    RmpADC_positive,
    rmpA2_positive,
    vf_factor_count_filtered
  ) %>%
  mutate(
    tree_order = match(
      SRR,
      tree$tip.label
    )
  ) %>%
  arrange(tree_order)

write_tsv(
  tree_meta,
  file.path(
    OUT,
    "01_Final1143_Tree_Metadata.tsv"
  )
)

# ============================================================
# 3. ELIGIBLE STs
# ============================================================

eligible <- read_tsv(
  ELIGIBLE_FILE,
  show_col_types=FALSE
)

eligible_STs <- eligible %>%
  distinct(ST)

cat(
  "Eligible STs for within-ST phylogenetic testing:",
  nrow(eligible_STs),
  "\n"
)

# ============================================================
# 4. WITHIN-ST PHYLOGENETIC STRUCTURE
#
# Compare patristic distance:
# same siderophore profile vs different profile
#
# Positive delta means:
# isolates sharing siderophore profile are closer.
# ============================================================

set.seed(20260812)

N_PERM <- 10000

test_ST <- function(st_name){

  z <- master %>%
    filter(ST == st_name) %>%
    semi_join(
      eligible %>% filter(ST == st_name),
      by=c("ST","siderophore_profile")
    )

  ids <- z$SRR

  if(length(ids) < 10){
    return(NULL)
  }

  st_tree <- keep.tip(
    tree,
    ids
  )

  D <- cophenetic.phylo(st_tree)

  meta <- z %>%
    select(
      SRR,
      siderophore_profile
    ) %>%
    distinct()

  prof <- setNames(
    meta$siderophore_profile,
    meta$SRR
  )

  labs <- prof[
    rownames(D)
  ]

  idx <- which(
    upper.tri(D),
    arr.ind=TRUE
  )

  distances <- D[idx]

  same <- labs[idx[,1]] ==
          labs[idx[,2]]

  if(
    sum(same) == 0 ||
    sum(!same) == 0
  ){
    return(NULL)
  }

  obs_same <- mean(
    distances[same]
  )

  obs_diff <- mean(
    distances[!same]
  )

  obs_delta <- obs_diff -
               obs_same

  perm_delta <- numeric(N_PERM)

  for(i in seq_len(N_PERM)){

    perm_labs <- sample(labs)

    same_perm <-
      perm_labs[idx[,1]] ==
      perm_labs[idx[,2]]

    if(
      sum(same_perm) == 0 ||
      sum(!same_perm) == 0
    ){
      perm_delta[i] <- NA_real_
    } else {

      perm_delta[i] <-
        mean(
          distances[!same_perm]
        ) -
        mean(
          distances[same_perm]
        )
    }
  }

  perm_delta <-
    perm_delta[
      !is.na(perm_delta)
    ]

  P_perm <- (
    sum(
      perm_delta >= obs_delta
    ) + 1
  ) / (
    length(perm_delta) + 1
  )

  tibble(
    ST=st_name,
    isolates=nrow(z),
    profiles=n_distinct(
      z$siderophore_profile
    ),
    mean_same_profile_distance=
      obs_same,
    mean_different_profile_distance=
      obs_diff,
    distance_delta=
      obs_delta,
    same_to_different_ratio=
      obs_same / obs_diff,
    permutation_P=P_perm
  )
}

phylo_tests <- lapply(
  eligible_STs$ST,
  test_ST
) %>%
  bind_rows() %>%
  mutate(
    FDR_P=p.adjust(
      permutation_P,
      method="BH"
    )
  ) %>%
  arrange(
    FDR_P,
    permutation_P
  )

write_tsv(
  phylo_tests,
  file.path(
    OUT,
    "02_WithinST_Siderophore_Phylogenetic_Structure.tsv"
  )
)

# ============================================================
# 5. PROFILE-LEVEL MRCA PURITY
#
# Measures how concentrated each ST-siderophore group
# is within its smallest containing clade.
# ============================================================

get_descendant_tips <- function(tree, node){

  Ntip <- length(
    tree$tip.label
  )

  children <- tree$edge[
    tree$edge[,1] == node,
    2
  ]

  descendants <- integer()

  recurse <- function(n){

    if(n <= Ntip){
      descendants <<-
        c(descendants,n)
    } else {
      kids <- tree$edge[
        tree$edge[,1] == n,
        2
      ]
      for(k in kids){
        recurse(k)
      }
    }
  }

  recurse(node)

  tree$tip.label[
    descendants
  ]
}

profile_results <- list()

counter <- 1

for(st_name in eligible_STs$ST){

  groups <- eligible %>%
    filter(ST == st_name)

  for(profile_name in groups$siderophore_profile){

    ids <- master %>%
      filter(
        ST == st_name,
        siderophore_profile ==
          profile_name
      ) %>%
      pull(SRR)

    if(length(ids) < 2)
      next

    node <- getMRCA(
      tree,
      ids
    )

    if(is.null(node))
      next

    descendants <-
      get_descendant_tips(
        tree,
        node
      )

    desc_meta <- master %>%
      filter(
        SRR %in% descendants
      )

    target_in_clade <-
      sum(
        desc_meta$ST == st_name &
        desc_meta$siderophore_profile ==
          profile_name
      )

    purity <-
      target_in_clade /
      nrow(desc_meta)

    ST_purity <-
      sum(
        desc_meta$ST ==
          st_name
      ) /
      nrow(desc_meta)

    profile_results[[counter]] <-
      tibble(
        ST=st_name,
        siderophore_profile=
          profile_name,
        group_isolates=
          length(ids),
        MRCA_clade_tips=
          nrow(desc_meta),
        target_isolates_in_MRCA=
          target_in_clade,
        profile_clade_purity=
          purity,
        ST_clade_purity=
          ST_purity
      )

    counter <- counter + 1
  }
}

profile_purity <- bind_rows(
  profile_results
) %>%
  arrange(
    ST,
    desc(profile_clade_purity)
  )

write_tsv(
  profile_purity,
  file.path(
    OUT,
    "03_ST_Siderophore_MRCA_Clade_Purity.tsv"
  )
)

# ============================================================
# 6. FOCAL LINEAGES
# ============================================================

focal <- c(
  "ST258",
  "ST147",
  "ST383",
  "ST16"
)

focal_results <- phylo_tests %>%
  filter(
    ST %in% focal
  )

write_tsv(
  focal_results,
  file.path(
    OUT,
    "04_Focal_Lineages_Phylogenetic_Structure.tsv"
  )
)

cat("\n========================================\n")
cat("WITHIN-ST PHYLOGENETIC RESULTS\n")
cat("========================================\n")

print(
  phylo_tests,
  n=Inf,
  width=Inf
)

cat("\n===== FOCAL LINEAGES =====\n")

print(
  focal_results,
  n=Inf,
  width=Inf
)

cat("\n===== FOCAL MRCA PURITY =====\n")

print(
  profile_purity %>%
    filter(
      ST %in% focal
    ),
  n=Inf,
  width=Inf
)

cat("\nFINAL 1143 PHYLOGENETIC STRUCTURE ANALYSIS COMPLETE\n")

