#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)

amr <- read_tsv(
  "data/02_AMR_Table.tsv",
  show_col_types = FALSE
) %>%
  mutate(
    KPC = str_detect(amr_gene_list, "blaKPC"),
    NDM = str_detect(amr_gene_list, "blaNDM"),
    OXA48_like = str_detect(amr_gene_list, "blaOXA-(48|181|232)"),
    VIM_IMP = str_detect(amr_gene_list, "bla(VIM|IMP)"),
    CTXM = str_detect(amr_gene_list, "blaCTX-M"),
    qnr = str_detect(amr_gene_list, "(^|;)qnr"),
    AAC6Ibcr = str_detect(amr_gene_list, "aac\\(6'\\)-Ib-cr"),
    gyrA_mut = str_detect(point_mutation_list, "(^|;)gyrA_"),
    parC_mut = str_detect(point_mutation_list, "(^|;)parC_"),
    mcr = str_detect(amr_gene_list, "(^|;)mcr[-0-9]"),
    colistin_regulatory = str_detect(
      point_mutation_list,
      "(^|;)(pmrA_|pmrB_|phoP_|phoQ_|mgrB_|crrB_)"
    )
  )

sid <- read_tsv(
  "data/04_Siderophore_Table.tsv",
  show_col_types = FALSE
) %>%
  select(SRR, siderophore_profile)

dat <- sid %>%
  inner_join(amr, by = "SRR")

major_profiles <- c(
  "none",
  "ybt",
  "ybt+iuc",
  "ybt+iuc+iro"
)

dat_major <- dat %>%
  filter(siderophore_profile %in% major_profiles)

markers <- c(
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

results <- lapply(markers, function(marker) {

  d <- dat_major %>%
    select(siderophore_profile, all_of(marker))

  tab <- table(
    factor(d$siderophore_profile, levels = major_profiles),
    d[[marker]]
  )

  test <- fisher.test(tab)

  out <- data.frame(
    Marker = marker,
    Fisher_P = test$p.value
  )

  for (p in major_profiles) {
    n_total <- sum(d$siderophore_profile == p)
    n_pos <- sum(
      d$siderophore_profile == p & d[[marker]],
      na.rm = TRUE
    )

    out[[paste0(p, "_positive")]] <- n_pos
    out[[paste0(p, "_total")]] <- n_total
    out[[paste0(p, "_percent")]] <- 100 * n_pos / n_total
  }

  out
})

results <- bind_rows(results) %>%
  mutate(
    FDR_BH = p.adjust(Fisher_P, method = "BH"),
    Significant_FDR_0_05 = FDR_BH < 0.05
  )

write_tsv(
  results,
  "results/siderophore/07_profile_AMR_marker_associations.tsv"
)

print(results)
