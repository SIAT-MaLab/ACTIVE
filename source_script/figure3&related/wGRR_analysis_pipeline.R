library(tidyverse)
edges_df <- read_csv("allpp_wGRR_Edges.csv", show_col_types = FALSE)
meta_df <- read_tsv("../../merged_phage_stats_taxonomy.tsv", show_col_types = FALSE)
meta_processed <- meta_df %>%
  mutate(Genome_ID = paste0(Contig_ID, "_", Start, "-", Stop)) %>%
  mutate(Activity_Score = as.numeric(Activity_Score),
         Activity_Score = replace_na(Activity_Score, 0),
         Status = if_else(Activity_Score >= 0.7, "Active", "Inactive")) %>%
  select(Genome_ID, Status, vFAM, vGENUS,
         host_Phylum, host_Class, host_Order, host_Family, host_Genus, host_Species, lifestyle) %>%
  distinct()

edges_annotated <- edges_df %>%
  left_join(meta_processed, by = c("query_genome" = "Genome_ID")) %>%
  rename_with(~paste0("Q_", .), c(Status, vFAM, vGENUS, host_Phylum, host_Class, host_Order, host_Family, host_Genus, host_Species, lifestyle)) %>%
  left_join(meta_processed, by = c("target_genome" = "Genome_ID")) %>%
  rename_with(~paste0("T_", .), c(Status, vFAM, vGENUS, host_Phylum, host_Class, host_Order, host_Family, host_Genus, host_Species, lifestyle)) %>%
  drop_na(Q_Status, T_Status) %>%
  filter(query_genome != target_genome) %>%
  mutate(
    Pair_Type = case_when(
      Q_Status == "Active" & T_Status == "Active" ~ "Active-Active",
      Q_Status == "Inactive" & T_Status == "Inactive" ~ "Inactive-Inactive",
      TRUE ~ "Active-Inactive"
    )
  )
my_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  )

step_size <- 0.01
ii_data <- edges_annotated %>% filter(Pair_Type == "Active-Active")

ii_freq_table <- ii_data %>%
  mutate(wGRR_Bin = cut(wGRR,
                        breaks = seq(0, 1, by = step_size),
                        include.lowest = TRUE,
                        right = FALSE)) %>%
  count(wGRR_Bin, name = "Pair_Count") %>%
  drop_na(wGRR_Bin) %>%
  mutate(Percentage = round(Pair_Count / sum(Pair_Count) * 100, 2)) %>%
  arrange(wGRR_Bin)

write_csv(ii_freq_table, sprintf("3_AA_wGRR_Frequency_Table_%.2f.csv", step_size))

p_ii_hist <- ggplot(ii_data, aes(x = wGRR)) +
  geom_histogram(binwidth = step_size, fill = "grey60", color = "black", alpha = 0.9, boundary = 0) +
  labs(title = sprintf("Active-Active wGRR Frequency (Binwidth = %s)", step_size),
       x = "wGRR Range", y = "Number of Pairs (Count)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.05)) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(sprintf("3_AA_wGRR_Histogram_%.2f.pdf", step_size), plot = p_ii_hist, width = 10, height = 5)
