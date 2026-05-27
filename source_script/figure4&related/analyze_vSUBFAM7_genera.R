

library(tidyverse)
library(ggstatsplot)
library(patchwork) 

INPUT_FILE <- "packaging_module_structure_vFAM7.csv"  
TARGET_SUBFAM <- "vSUBFAM_7"
MIN_GENOME_PER_GENUS <- 3 
df <- read.csv(INPUT_FILE)

df_clean <- df %>%
  filter(vSUBFAM == TARGET_SUBFAM) %>%
  filter(!is.na(Gap_HNH_TerS_Count) & !is.na(Gap_TerS_TerL_Count)) %>%
  filter(!is.na(vGENUS) & vGENUS != "") %>%
  group_by(vGENUS) %>%
  filter(n() >= MIN_GENOME_PER_GENUS) %>%
  ungroup()

genus_count <- length(unique(df_clean$vGENUS))
cat("------------------------------------------------\n")
cat("Target Subfamily:", TARGET_SUBFAM, "\n")
cat("Genera included (n >=", MIN_GENOME_PER_GENUS, "):", genus_count, "\n")
cat("Total genomes analyzed:", nrow(df_clean), "\n")
cat("------------------------------------------------\n")

if (nrow(df_clean) == 0) {
  stop("Error: No data found after filtering. Check if vSUBFAM_7 exists in your CSV.")
}

print("Generating Plot 1: AB Gap Analysis (by vGENUS)...")

p1 <- ggbetweenstats(
  data = df_clean,
  x = vGENUS,
  y = Gap_HNH_TerS_Count,
  type = "nonparametric", # Kruskal-Wallis
  title = paste0("HNH-TerS Gap in ", TARGET_SUBFAM, " (by Genus)"),
  xlab = "Viral Genus",
  ylab = "Inserted CDS Count (A to B)",
  package = "ggthemes",
  palette = "Tableau_10",
  plot.type = "box",
  results.subtitle = TRUE,
  centrality.plotting = TRUE,
  pairwise.display = "none",
  p.adjust.method = "fdr"
) +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  plot.title = element_text(size = 14, face = "bold")
)

print("Generating Plot 2: BC Gap Analysis (by vGENUS)...")

p2 <- ggbetweenstats(
  data = df_clean,
  x = vGENUS,
  y = Gap_TerS_TerL_Count,
  type = "nonparametric",
  title = paste0("TerS-TerL Gap in ", TARGET_SUBFAM, " (by Genus)"),
  xlab = "Viral Genus",
  ylab = "Inserted CDS Count (B to C)",
  package = "ggthemes",
  palette = "Classic_10",
  plot.type = "box",
  results.subtitle = TRUE,
  centrality.plotting = TRUE,
  pairwise.display = "none",
  p.adjust.method = "fdr"
) +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  plot.title = element_text(size = 14, face = "bold")
)

ggsave(paste0(TARGET_SUBFAM, "_Genus_Stats_AB.pdf"), p1, width = 12, height = 8)
ggsave(paste0(TARGET_SUBFAM, "_Genus_Stats_BC.pdf"), p2, width = 12, height = 8)

print("Done! Files generated:")
print(paste0(TARGET_SUBFAM, "_Genus_Stats_AB.pdf"))
print(paste0(TARGET_SUBFAM, "_Genus_Stats_BC.pdf"))

