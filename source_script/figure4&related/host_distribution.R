
library(tidyverse)
library(ggplot2)
library(scales)      
library(patchwork)   #
PHYLUM_COLORS <- c(
  "Bacillota_A"       = "#dce8ef",
  "Bacillota_C"       = "#add2e5",
  "Bacillota_I"       = "#34666b",
  "Bacillota"         = "#7abbce",
  "Bacillota_B"       = "#769499",
  "Actinomycetota"    = "#f5cfa6",
  "Bacteroidota"      = "#efe3ef",
  "Verrucomicrobiota" = "#d4d68a",
  "Pseudomonadota"    = "#e68f9f",
  "Desulfobacterota"  = "#d5d5d6",
  "Fusobacteriota"    = "#b7a3c9"
)


PIE_BORDER_COLOR <- "black"
PIE_BORDER_WIDTH <- 0.4

input_file <- "../../../merged_phage_stats_taxonomy.tsv"
meta_data <- read_tsv(input_file, show_col_types = FALSE)

vfam7_data <- meta_data %>% filter(vFAM == "vFAM_7")

pie_phylum <- vfam7_data %>%
  filter(!is.na(Phylum) & Phylum != "" & Phylum != "None") %>%
  count(Phylum) %>%
  arrange(desc(n)) %>%
  mutate(
    prop = n / sum(n),
    label_text = paste0(Phylum, " (n=", n, ", ", percent(prop, accuracy = 0.1), ")")
  )

pie_phylum$Phylum <- factor(pie_phylum$Phylum, levels = pie_phylum$Phylum)

missing_phyla <- setdiff(pie_phylum$Phylum, names(PHYLUM_COLORS))
if(length(missing_phyla) > 0) {
  extra_colors <- setNames(rep("#95a5a6", length(missing_phyla)), missing_phyla)
  PHYLUM_COLORS <- c(PHYLUM_COLORS, extra_colors)
}


pie_family <- vfam7_data %>%
  filter(!is.na(Family) & Family != "" & Family != "None") %>%
  count(Family) %>%
  mutate(Family_Group = ifelse(n >= 10, Family, "Other")) %>%
  group_by(Family_Group) %>%
  summarise(n = sum(n), .groups = 'drop') %>%
  arrange(desc(n)) %>%
  mutate(
    prop = n / sum(n),
    label_text = paste0(Family_Group, " (n=", n, ", ", percent(prop, accuracy = 0.1), ")")
  )

family_levels <- setdiff(pie_family$Family_Group, "Other")
if ("Other" %in% pie_family$Family_Group) {
  family_levels <- c(family_levels, "Other")
}
pie_family$Family_Group <- factor(pie_family$Family_Group, levels = family_levels)
pie_family <- pie_family %>% arrange(Family_Group)

num_main_families <- length(setdiff(family_levels, "Other"))
family_colors <- hue_pal()(num_main_families)
names(family_colors) <- setdiff(family_levels, "Other")
if ("Other" %in% family_levels) {
  family_colors["Other"] <- "#bdc3c7" 
}


p1 <- ggplot(pie_phylum, aes(x = "", y = n, fill = Phylum)) +
  geom_bar(stat = "identity", width = 1, color = PIE_BORDER_COLOR, linewidth = PIE_BORDER_WIDTH) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(
    values = PHYLUM_COLORS,
    breaks = pie_phylum$Phylum, 
    labels = pie_phylum$label_text
  ) +
  theme_void() +
  labs(title = "A. Host Phylum Distribution", fill = "Phylum") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )

p2 <- ggplot(pie_family, aes(x = "", y = n, fill = Family_Group)) +
  geom_bar(stat = "identity", width = 1, color = PIE_BORDER_COLOR, linewidth = PIE_BORDER_WIDTH) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(
    values = family_colors,
    breaks = pie_family$Family_Group, 
    labels = pie_family$label_text
  ) +
  theme_void() +
  labs(title = "B. Host Family Distribution (Threshold >= 10)", fill = "Family") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )

final_plot <- p1 + p2 + 
  plot_annotation(
    title = "Host Taxonomy Distribution for vFAM_7 Phages",
    subtitle = paste0("Total occurrences analyzed: n = ", nrow(vfam7_data)),
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 20))
    )
  )

output_file <- "vFAM7_Host_Distribution_Combined.pdf"
ggsave(output_file, final_plot, width = 16, height = 7, bg = "white", device = cairo_pdf)

