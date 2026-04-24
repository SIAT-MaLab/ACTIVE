
library(tidyverse)
library(ggstatsplot)

df <- read.csv("packaging_module_structure_vFAM7.csv")
df_clean <- df %>%
  filter(!is.na(vSUBFAM) & vSUBFAM != "") %>%
  filter(!is.na(Gap_HNH_TerS_Count) & !is.na(Gap_TerS_TerL_Count)) %>%
  group_by(vSUBFAM) %>%
  filter(n() >= 3) %>%
  ungroup()
df_plot <- df_clean

p1 <- ggbetweenstats(
  data = df_plot,
  x = vSUBFAM,
  y = Gap_HNH_TerS_Count,
  type = "nonparametric", 
  title = "Number of CDS Inserted between HNH and TerS (AB)",
  xlab = "Viral Subfamily",
  ylab = "Inserted CDS Count (A to B)",
  package = "ggthemes",
  palette = "Tableau_10",
  plot.type = "boxviolin",
  results.subtitle = TRUE, 
  centrality.plotting = TRUE, 
  p.adjust.method = "fdr"
) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Stats_AB_Gap_vSUBFAM.pdf", p1, width = 14, height = 8)

p2 <- ggbetweenstats(
  data = df_plot,
  x = vSUBFAM,
  y = Gap_TerS_TerL_Count,
  type = "nonparametric",
  title = "Number of CDS Inserted between TerS and TerL (BC)",
  xlab = "Viral Subfamily",
  ylab = "Inserted CDS Count (B to C)",
  package = "ggthemes",
  palette = "Classic_10",
  plot.type = "boxviolin",
  results.subtitle = TRUE,
  centrality.plotting = TRUE,
  p.adjust.method = "fdr"
) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Stats_BC_Gap_vSUBFAM.pdf", p2, width = 14, height = 8)

print("Analysis complete. Check Stats_AB_Gap_vSUBFAM.pdf and Stats_BC_Gap_vSUBFAM.pdf")

