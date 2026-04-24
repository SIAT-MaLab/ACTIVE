library(tidyverse)
library(ggpubr)

tax_file <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/merged_phage_stats_taxonomy.tsv"
map_file <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/genome_to_vOTU.tsv"

tax_df <- read_tsv(tax_file, show_col_types = FALSE)

vfam10_status <- tax_df %>%
  filter(vFAM == "vFAM_10") %>%
  group_by(vOTU) %>%
  summarise(
    Is_Active = any(!is.na(Activity_Score) & Activity_Score >= 0.7),
    .groups = "drop"
  ) %>%
  mutate(Status = ifelse(Is_Active, "Active", "Inactive"))

map_df <- read_tsv(map_file, col_names = c("Genome", "vOTU", "Rep_Genome"), show_col_types = FALSE)

vfam10_reps <- map_df %>%
  filter(vOTU %in% vfam10_status$vOTU) %>%
  select(vOTU, Rep_Genome) %>%
  distinct() %>%
  left_join(vfam10_status, by = "vOTU") %>%
  mutate(coord_string = str_extract(Rep_Genome, "\\d+-\\d+$")) %>%
  separate(coord_string, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  mutate(Length_bp = end - start + 1) %>%
  filter(!is.na(Length_bp))

plot_data <- bind_rows(
  vfam10_reps %>% mutate(Group = Status),        
  vfam10_reps %>% mutate(Group = "All")          
) %>%
  mutate(Group = factor(Group, levels = c("Active", "Inactive", "All")))

medians <- plot_data %>%
  group_by(Group) %>%
  summarise(
    Median_Length = median(Length_bp),
    Count = n(),
    .groups = "drop"
  )


print(as.data.frame(medians))

p <- ggplot(plot_data, aes(x = Group, y = Length_bp, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("Active" = "#e31a1c", "Inactive" = "#1f78b4", "All" = "#33a02c")) +
  labs(
    title = "Comparison of Representative Genome Lengths in vFAM_10",
    subtitle = "Active defined as vOTU with max Activity_Score >= 0.7",
    x = "vOTU Activity Status",
    y = "Genome Length (bp)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(comparisons = list(c("Active", "Inactive")), method = "wilcox.test", label = "p.signif")

out_img <- "vFAM10_Length_Comparison.pdf"
ggsave(out_img, p, width = 7, height = 6)
