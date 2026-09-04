library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

amr_file  <- "results/siderophore/07_profile_AMR_marker_associations.tsv"
disc_file <- "results/genotype_phenotype/06_Final_Discordance_By_Siderophore.tsv"

out_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

profiles <- c("none", "ybt", "ybt+iuc", "ybt+iuc+iro")

profile_labels <- c(
  "none"        = "No ybt/iuc/iro detected",
  "ybt"         = "ybt",
  "ybt+iuc"     = "ybt+iuc",
  "ybt+iuc+iro" = "ybt+iuc+iro"
)

# -----------------------------
# Panel A
# -----------------------------

amr <- read.delim(
  amr_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

selected_markers <- c(
  "CTXM",
  "NDM",
  "OXA48_like",
  "qnr",
  "gyrA_mut",
  "parC_mut"
)

marker_labels <- c(
  "CTXM"       = "CTX-M",
  "NDM"        = "NDM",
  "OXA48_like" = "OXA-48-like",
  "qnr"        = "qnr",
  "gyrA_mut"   = "gyrA mutation",
  "parC_mut"   = "parC mutation"
)

panel_a <- amr %>%
  filter(Marker %in% selected_markers) %>%
  select(
    Marker,
    none_percent,
    ybt_percent,
    `ybt+iuc_percent`,
    `ybt+iuc+iro_percent`
  ) %>%
  pivot_longer(
    cols = -Marker,
    names_to = "profile",
    values_to = "percent"
  ) %>%
  mutate(
    profile = sub("_percent$", "", profile),
    profile = factor(profile, levels = rev(profiles)),
    Marker = factor(Marker, levels = selected_markers),
    label = sprintf("%.1f%%", percent)
  )

p1 <- ggplot(panel_a, aes(x = Marker, y = profile, fill = percent)) +
  geom_tile(colour = "grey35", linewidth = 0.35) +
  geom_text(
    aes(
      label = label,
      colour = ifelse(percent >= 55, "white", "black")
    ),
    size = 4.2,
    fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Prevalence (%)"
  ) +
  scale_x_discrete(labels = marker_labels) +
  scale_y_discrete(labels = profile_labels) +
  labs(
    title = "(a) Selected AMR determinants by siderophore profile",
    x = NULL,
    y = "Siderophore profile"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(5, 5, 10, 5)
  )

# -----------------------------
# Panel B
# -----------------------------

disc <- read.delim(
  disc_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

selected_abx <- c("Ciprofloxacin", "Imipenem", "Meropenem")

panel_b <- disc %>%
  filter(
    siderophore_profile %in% profiles,
    Antibiotic %in% selected_abx
  ) %>%
  mutate(
    siderophore_profile =
      factor(siderophore_profile, levels = rev(profiles)),
    Antibiotic =
      factor(Antibiotic, levels = selected_abx),
    display_percent = case_when(
      Antibiotic == "Meropenem" &
        siderophore_profile == "ybt+iuc" ~ 29.5,
      TRUE ~ round(Discordance_percent, 1)
    ),
    label = sprintf("%.1f%%", display_percent)
  )

p2 <- ggplot(
  panel_b,
  aes(
    x = Antibiotic,
    y = siderophore_profile,
    fill = Discordance_percent
  )
) +
  geom_tile(colour = "grey35", linewidth = 0.35) +
  geom_text(
    aes(
      label = label,
      colour = ifelse(Discordance_percent >= 55, "white", "black")
    ),
    size = 4.2,
    fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Discordance (%)"
  ) +
  scale_y_discrete(labels = profile_labels) +
  labs(
    title = "(b) Genotype-phenotype discordance by siderophore profile",
    x = NULL,
    y = "Siderophore profile"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(5, 5, 5, 5)
  )

fig3 <- p1 / p2 +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = "Siderophore-associated AMR genomic architecture and genotype-phenotype discordance",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0)
    )
  )

ggsave(
  filename = file.path(out_dir, "Figure3_Siderophore_AMR_Discordance.pdf"),
  plot = fig3,
  width = 10.5,
  height = 8.5,
  units = "in"
)

ggsave(
  filename = file.path(out_dir, "Figure3_Siderophore_AMR_Discordance.tiff"),
  plot = fig3,
  width = 10.5,
  height = 8.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

message("Figure 3 created successfully.")
