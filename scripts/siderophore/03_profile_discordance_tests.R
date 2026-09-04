#!/usr/bin/env Rscript

input_file <- "results/genotype_phenotype/06_Final_Discordance_By_Siderophore.tsv"
output_file <- "results/siderophore/06_profile_discordance_fisher.tsv"

x <- read.delim(input_file, stringsAsFactors = FALSE)

major_profiles <- c(
  "none",
  "ybt",
  "ybt+iuc",
  "ybt+iuc+iro"
)

x <- x[x$siderophore_profile %in% major_profiles, ]
x$Concordant <- x$Evaluated - x$Discordant

antibiotics <- sort(unique(x$Antibiotic))

results <- lapply(antibiotics, function(ab) {

  d <- x[x$Antibiotic == ab, ]

  d <- d[match(major_profiles, d$siderophore_profile), ]

  m <- cbind(
    Discordant = d$Discordant,
    Concordant = d$Concordant
  )

  rownames(m) <- d$siderophore_profile

  test <- fisher.test(m)

  data.frame(
    Antibiotic = ab,
    Evaluated = sum(d$Evaluated),
    Discordant = sum(d$Discordant),
    Fisher_P = test$p.value
  )
})

results <- do.call(rbind, results)

results$FDR_BH <- p.adjust(
  results$Fisher_P,
  method = "BH"
)

results$Significant_FDR_0_05 <- results$FDR_BH < 0.05

write.table(
  results,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(results)
