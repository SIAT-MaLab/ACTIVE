
library(tidyverse)
library(ggsci)

df_effects <- read_tsv("summary_all_taxa_effects.tsv", show_col_types = FALSE)

plot_data_cross <- df_effects %>%
  filter(Trait == "Induction Rate") %>%
  filter(Effective_N >= 3) %>%
  mutate(Display_Label = paste0(Short_Name, " (", Level, ")")) %>%
  arrange(Effect)
n_rows <- nrow(plot_data_cross)

if (n_rows > 30) {
  plot_ready <- plot_data_cross %>% slice(c(1:15, (n()-14):n()))
} else {
  plot_ready <- plot_data_cross
}

level_order <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
plot_ready$Level <- factor(plot_ready$Level, levels = level_order)

p_cross <- ggplot(plot_ready, aes(x = reorder(Display_Label, Effect), y = Effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_segment(aes(x = Display_Label, xend = Display_Label, y = 0, yend = Effect), color = "grey80") +
  geom_point(aes(size = Effective_N, color = Level), alpha = 0.8) +
  coord_flip() +

  scale_size_continuous(
    range = c(1, 15),  
    breaks = c(10, 100, 1000, 5000) 
  ) +

  scale_color_npg() +
  labs(title = "Top & Bottom Taxa for Induction Rate (Cross-Level)",
       subtitle = "Bubbles sized by Sample Size (Range 1-15)",
       x = "", 
       y = "Effect Size (Log-Odds: >0 Active, <0 Silent)",
       size = "Sample Size",
       color = "Taxonomic Level") +
  
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10), 
    legend.position = "right",
    legend.box = "vertical"
  )

ggsave("plot_induction_cross_level_bubble_v2.pdf", p_cross, width = 11, height = 12)

