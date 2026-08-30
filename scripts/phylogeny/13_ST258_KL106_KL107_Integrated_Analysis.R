

library(readr)
library(dplyr)
library(tibble)

# Repository root
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), "..", ".."))

KFILE <- file.path(
  base_dir, "data", "ST258_Klocus_FINAL1143.tsv"
)

GPFILE <- file.path(
  base_dir, "data",
  "03_Isolate_Level_Siderophore_AMR_Genotype_Phenotype.tsv"
)

output_dir <- file.path(base_dir, "results", "phylogeny")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

OUTFILE <- file.path(
  output_dir,
  "07_ST258_KL106_KL107_Integrated_Analysis.tsv"
)

# ------------------------------------------------------------
# ST258 KL106 / KL107 cohort
# ------------------------------------------------------------

k <- read_tsv(
  KFILE,
  show_col_types = FALSE
) %>%
  filter(K_locus %in% c("KL106", "KL107"))

# ------------------------------------------------------------
# Helper function for 2x2 Fisher test
# ------------------------------------------------------------

run_test <- function(
  feature,
  a_pos, a_neg,
  b_pos, b_neg,
  a_label = "KL106",
  b_label = "KL107"
) {

  x <- matrix(
    c(
      a_pos, a_neg,
      b_pos, b_neg
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      K_locus = c(a_label, b_label),
      Status = c("Positive_or_R", "Negative_or_S")
    )
  )

  f <- fisher.test(x)

  tibble(
    Feature = feature,
    KL106_positive_or_R = a_pos,
    KL106_negative_or_S = a_neg,
    KL107_positive_or_R = b_pos,
    KL107_negative_or_S = b_neg,
    KL106_percent = 100 * a_pos / (a_pos + a_neg),
    KL107_percent = 100 * b_pos / (b_pos + b_neg),
    OR_KL106_vs_KL107 = unname(f$estimate),
    CI_low = f$conf.int[1],
    CI_high = f$conf.int[2],
    Fisher_P = f$p.value
  )
}

# ------------------------------------------------------------
# YERSINIABACTIN
# ------------------------------------------------------------

ybt_106_pos <- sum(k$K_locus == "KL106" & k$Yersiniabactin != "-")
ybt_106_neg <- sum(k$K_locus == "KL106" & k$Yersiniabactin == "-")

ybt_107_pos <- sum(k$K_locus == "KL107" & k$Yersiniabactin != "-")
ybt_107_neg <- sum(k$K_locus == "KL107" & k$Yersiniabactin == "-")

ybt_res <- run_test(
  "Yersiniabactin",
  ybt_106_pos,
  ybt_106_neg,
  ybt_107_pos,
  ybt_107_neg
)

# ------------------------------------------------------------
# COLIBACTIN
# ------------------------------------------------------------

clb_106_pos <- sum(k$K_locus == "KL106" & k$Colibactin != "-")
clb_106_neg <- sum(k$K_locus == "KL106" & k$Colibactin == "-")

clb_107_pos <- sum(k$K_locus == "KL107" & k$Colibactin != "-")
clb_107_neg <- sum(k$K_locus == "KL107" & k$Colibactin == "-")

clb_res <- run_test(
  "Colibactin",
  clb_106_pos,
  clb_106_neg,
  clb_107_pos,
  clb_107_neg
)

# ------------------------------------------------------------
# PHENOTYPIC RESISTANCE
# ------------------------------------------------------------

gp <- read_tsv(
  GPFILE,
  show_col_types = FALSE
) %>%
  filter(Phenotype %in% c("R", "S")) %>%
  distinct(SRR, Antibiotic, .keep_all = TRUE) %>%
  inner_join(
    k %>% select(SRR, K_locus),
    by = "SRR"
  )

get_pheno_counts <- function(drug) {

  z <- gp %>%
    filter(Antibiotic == drug)

  a_pos <- sum(z$K_locus == "KL106" & z$Phenotype == "R")
  a_neg <- sum(z$K_locus == "KL106" & z$Phenotype == "S")

  b_pos <- sum(z$K_locus == "KL107" & z$Phenotype == "R")
  b_neg <- sum(z$K_locus == "KL107" & z$Phenotype == "S")

  run_test(
    paste0(drug, "_phenotypic_resistance"),
    a_pos,
    a_neg,
    b_pos,
    b_neg
  )
}

cip_res <- get_pheno_counts("Ciprofloxacin")
ipm_res <- get_pheno_counts("Imipenem")
mem_res <- get_pheno_counts("Meropenem")

# ------------------------------------------------------------
# COMBINE + FDR
# ------------------------------------------------------------

final_res <- bind_rows(
  ybt_res,
  clb_res,
  cip_res,
  ipm_res,
  mem_res
) %>%
  mutate(
    BH_FDR = p.adjust(
      Fisher_P,
      method = "BH"
    )
  )

write_tsv(
  final_res,
  OUTFILE
)

cat("\n============================================\n")
cat("ST258 KL106 vs KL107 INTEGRATED ANALYSIS\n")
cat("============================================\n\n")

print(
  final_res,
  n = Inf,
  width = Inf
)

cat("\nSaved to:\n", OUTFILE, "\n")

