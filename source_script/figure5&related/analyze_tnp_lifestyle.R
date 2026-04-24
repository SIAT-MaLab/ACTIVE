
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

BASE_DIR <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats"
META_FILE <- file.path(BASE_DIR, "allpp_checkv_quality/hqphage_metadata.tsv")
LIFESTYLE_FILE <- file.path(BASE_DIR, "allpp_annot/pfamscan/Phage_Lifestyle.tsv")

OUTPUT_TSV <- "Caudo_Mu_Lifestyle_Stats.tsv"
OUTPUT_PLOT <- "Caudo_Mu_Lifestyle_Distribution.pdf"

df_meta <- read_tsv(META_FILE, show_col_types = FALSE)

caudo_phages <- df_meta %>%
  filter(str_detect(ictv_Class, "(?i)Caudo")) %>%
  pull(Phage_Name) %>%
  unique()

caudo_phages_clean <- str_replace_all(caudo_phages, "-", "_")

df_life <- read_tsv(LIFESTYLE_FILE, show_col_types = FALSE)

df_life <- df_life %>%
  mutate(contig_clean = str_replace(Contig, "_(\\d+)-(\\d+)$", "_\\1_\\2"))

df_life_caudo <- df_life %>%
  filter(Contig %in% caudo_phages | contig_clean %in% caudo_phages_clean)

df_mu_stats <- df_life_caudo %>%
  filter(Analyzed_Mu_Lifestyle != "Unknown" & !is.na(Analyzed_Mu_Lifestyle)) %>%
  group_by(Analyzed_Mu_Lifestyle) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(desc(Count))

if (nrow(df_mu_stats) == 0) {
  quit(save = "no")
}
write_tsv(df_mu_stats, OUTPUT_TSV)


total_analyzed <- sum(df_mu_stats$Count)

p <- ggplot(df_mu_stats, aes(x = reorder(Analyzed_Mu_Lifestyle, Percentage), y = Percentage, fill = Analyzed_Mu_Lifestyle)) +
  geom_col(color = "black", alpha = 0.8, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%% (%d)", Percentage, Count)), 
            hjust = -0.1, size = 4.5, fontface = "bold", color = "black") +
  coord_flip() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_fill_brewer(palette = "Set2") + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Distribution of Analyzed Mu Lifestyles in Caudoviricetes",
    subtitle = sprintf("Excluding 'Unknown' categories (Total valid contigs: %d)", total_analyzed),
    x = "Mu-like Lifestyle Category",
    y = "Percentage (%)"
  )

ggsave(OUTPUT_PLOT, plot = p, width = 9, height = 6, dpi = 300)

print(df_mu_stats %>% mutate(Percentage = sprintf("%.2f%%", Percentage)))
cat("=======================================================\n")
