
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

BASE_DIR <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats"
ANNOT_FILE <- file.path(BASE_DIR, "allpp_annot/allpharokka_annot.tsv")
META_FILE <- file.path(BASE_DIR, "allpp_checkv_quality/hqphage_metadata.tsv")

OUTPUT_TSV <- "Transposase_Adjacent_Operon_Genes_PerContig.tsv"
OUTPUT_PLOT <- "Transposase_Adjacent_Top10.pdf" 

df_meta <- read_tsv(META_FILE, show_col_types = FALSE)

caudo_phages <- df_meta %>%
  filter(str_detect(ictv_Class, "(?i)Caudo")) %>%
  pull(Phage_Name) %>%
  unique()

df_cds <- read_tsv(ANNOT_FILE, show_col_types = FALSE)
df_cds <- df_cds %>%
  mutate(contig_clean = str_replace(contig, "_(\\d+)-(\\d+)$", "_\\1_\\2"))

caudo_phages_clean <- str_replace_all(caudo_phages, "-", "_")

df_cds_filtered <- df_cds %>%
  filter(contig %in% caudo_phages | contig_clean %in% caudo_phages_clean)

df_processed <- df_cds_filtered %>%
  mutate(
    phys_start = pmin(start, stop),
    annot_lower = tolower(annot),
    is_transposase = str_detect(annot_lower, "transposase|transposon")
  ) %>%
  arrange(contig, phys_start)

df_adj <- df_processed %>%
  group_by(contig) %>%
  mutate(
    prev_is_tnp = lag(is_transposase, default = FALSE),
    prev_frame  = lag(frame, default = "NA"),
    next_is_tnp = lead(is_transposase, default = FALSE),
    next_frame  = lead(frame, default = "NA")
  ) %>%
  filter(
    is_transposase == FALSE & (
      (prev_is_tnp == TRUE & prev_frame == frame) |
      (next_is_tnp == TRUE & next_frame == frame)
    )
  ) %>%
  ungroup() %>%
  distinct(gene, .keep_all = TRUE) 

if (nrow(df_adj) == 0) {
  quit(save = "no")
}

df_unique_per_contig <- df_adj %>%
  distinct(contig, phrog, category, annot)

total_valid_contigs <- n_distinct(df_adj$contig)

df_stats <- df_unique_per_contig %>%
  group_by(phrog, category, annot) %>%
  summarise(Contig_Count = n(), .groups = "drop") %>%
  mutate(Genome_Frequency_Pct = (Contig_Count / total_valid_contigs) * 100) %>%
  arrange(desc(Contig_Count))

write_tsv(df_stats, OUTPUT_TSV)

top10 <- df_stats %>% head(10)

top10 <- top10 %>%
  mutate(
    short_annot = str_trunc(annot, 40, "right"),
    plot_label = sprintf("%s\n(%s)", phrog, short_annot)
  )

p <- ggplot(top10, aes(x = reorder(plot_label, Genome_Frequency_Pct), y = Genome_Frequency_Pct, fill = category)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  coord_flip() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Top 10 Co-directional Genes Adjacent to Transposase",
    subtitle = sprintf("Calculated per genome across %d valid Caudoviricetes contigs", total_valid_contigs),
    x = "PHROG & Annotation",
    y = "Genome Presence Frequency (%)",
    fill = "Phrog Category"
  )

ggsave(OUTPUT_PLOT, plot = p, width = 10, height = 7, dpi = 300)

print(top10 %>% select(phrog, category, annot, Contig_Count, Genome_Frequency_Pct))
