
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork) 
})

META_FILE    <- "hqphage_metadata.tsv"
PC_MAP_FILE  <- "allpp_PC_map.tsv"
ANNOT_FILE   <- "mmseqs_repseq_annot/final/annotations-Protein_mmseqs_out_rep_seq/final_annotation_summary.tsv"
RESULT_FILE  <- "Final_TwoStep_BestHit.csv"

OUT_PREFIX   <- "Plot_Top10_BestHit"

meta    <- fread(META_FILE)
pc_map  <- fread(PC_MAP_FILE)
annot   <- fread(ANNOT_FILE)
res_df  <- fread(RESULT_FILE)

meta_clean <- meta %>%
  filter(grepl("Caudoviricetes", ictv_Class)) %>%
  mutate(is_active = if_else(!is.na(Activity_Score) & Activity_Score >= 0.7, 1, 0))

otu_stats <- meta_clean %>%
  group_by(vOTU) %>%
  summarise(
    total_genomes = n(),
    active_count  = sum(is_active),
    activity_rate = sum(is_active) / n(), 
    .groups = "drop"
  ) %>%
  filter(total_genomes >= 5)

valid_votus <- otu_stats$vOTU

top10_features <- res_df %>%
  filter(GLMM_FDR < 0.05, GLMM_Est > 0) %>%
  arrange(desc(GLMM_Est)) %>%
  head(10) %>%
  mutate(
    Clean_Product = str_trunc(Product, 50, "right"),
    Label = paste0(Clean_Product, " (", Feature, ")")
  )

top10_features$Label <- factor(top10_features$Label, levels = rev(top10_features$Label))
target_feats <- top10_features$Feature

feat_map <- pc_map %>%
  mutate(Phage_Name = str_remove(protein_id, "_CDS_[0-9]+$")) %>%
  inner_join(meta_clean %>% select(Phage_Name, vOTU), by = "Phage_Name") %>%
  filter(vOTU %in% valid_votus) %>%
  inner_join(annot, by = c("rep_seq_id" = "target")) %>%
  filter(best_hit %in% target_feats) %>%
  distinct(vOTU, Feature = best_hit) %>%
  mutate(has_feature = TRUE)

plot_summary_list <- list() # 用于画图的汇总数据
detail_raw_list   <- list() # 用于导出的明细数据

for(f in target_feats) {
  curr_prod <- top10_features %>% filter(Feature == f) %>% pull(Product)
  votus_with <- feat_map %>% filter(Feature == f) %>% pull(vOTU)
  df_detail <- otu_stats %>%
    mutate(
      Feature = f,
      Product = curr_prod,
      Group = if_else(vOTU %in% votus_with, "With", "Without")
    ) %>%
    select(Feature, Product, vOTU, Group, activity_rate, total_genomes, active_count)
  
  detail_raw_list[[f]] <- df_detail
  
  stats_with <- df_detail %>% filter(Group == "With") %>%
    summarise(
      mean_rate = sum(active_count) / sum(total_genomes), 
      n = n()
    ) %>% mutate(group = "With")
  
  stats_without <- df_detail %>% filter(Group == "Without") %>%
    summarise(
      mean_rate = sum(active_count) / sum(total_genomes), 
      n = n()
    ) %>% mutate(group = "Without")
  
  plot_summary_list[[f]] <- bind_rows(stats_with, stats_without) %>% mutate(Feature = f)
}

all_details <- bind_rows(detail_raw_list)
all_details <- all_details %>%
  arrange(factor(Feature, levels = target_feats), desc(Group), desc(activity_rate))
detail_file <- paste0(OUT_PREFIX, "_Detail_Statistics.csv")
write.csv(all_details, detail_file, row.names = FALSE)

plot_df_wide <- bind_rows(plot_summary_list) %>%
  left_join(top10_features %>% select(Feature, Label), by = "Feature") %>%
  pivot_wider(names_from = group, values_from = c(mean_rate, n))

p_left <- ggplot(plot_df_wide, aes(y = Label)) +
  geom_segment(aes(x = mean_rate_Without, xend = mean_rate_With, yend = Label),
               color = "grey80", linewidth = 1.5) +
  geom_point(aes(x = mean_rate_Without), size = 6, color = "grey60") +
  geom_point(aes(x = mean_rate_With), size = 8, color = "#377eb8") +
  scale_x_continuous(labels = percent_format(accuracy = 1), 
                     name = "Mean Activity Rate of vOTUs") +
  labs(title = "Observed Activity", subtitle = "Grey: Without Feature | Blue: With Feature") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(r = 20)
  )

p_right <- ggplot(top10_features, aes(y = Label, x = GLMM_Est)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 8, color = "#e41a1c") +
  scale_x_continuous(name = "GLMM Estimate (Log-Odds)", limits = c(1, NA)) +
  labs(title = "Model Effect Size", subtitle = "Positive = Pro-Active") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(l = 5)
  )

combined_plot <- p_left + p_right + plot_layout(widths = c(1.3, 1))

out_png <- paste0(OUT_PREFIX, ".png")
out_pdf <- paste0(OUT_PREFIX, ".pdf")

ggsave(out_pdf, combined_plot, width = 15, height = 7)
ggsave(out_png, combined_plot, width = 15, height = 7, dpi = 300)

