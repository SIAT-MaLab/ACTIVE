library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(scales)
library(patchwork)

META_FILE    <- "hqphage_metadata.tsv"
PC_MAP_FILE  <- "allpp_PC_map.tsv"
ANNOT_FILE   <- "mmseqs_repseq_annot/final/annotations-Protein_mmseqs_out_rep_seq/final_annotation_summary.tsv"
OUT_PREFIX   <- "Final_Analysis_BestHit_"

meta    <- fread(META_FILE)
pc_map  <- fread(PC_MAP_FILE)
annot   <- fread(ANNOT_FILE)

get_mode <- function(x) {
  ux <- unique(na.omit(x))
  if(length(ux) == 0) return("Unknown")
  ux[which.max(tabulate(match(x, ux)))]
}

meta_clean <- meta %>%
  filter(grepl("Caudoviricetes", ictv_Class)) %>%
  mutate(is_active = if_else(!is.na(Activity_Score) & Activity_Score >= 0.7, 1, 0))

otu_stats <- meta_clean %>%
  group_by(vOTU) %>%
  summarise(
    total_genomes = n(),
    active_count = sum(is_active),
    inactive_count = n() - sum(is_active),
    activity_rate = sum(is_active) / n(),
    host_phylum = get_mode(Phylum),
    host_class  = get_mode(Class),
    vFAM        = get_mode(vFAM),
    vGENUS      = get_mode(vGENUS),
    location_type = get_mode(Phage_Location_Type),
    .groups = "drop"
  ) %>%
  filter(total_genomes >= 5)

target_taxa_col <- "host_phylum"
top_taxa <- names(sort(table(otu_stats[[target_taxa_col]]), decreasing = TRUE))[1:15]

otu_stats <- otu_stats %>%
  mutate(host_taxa_reg = if_else(.data[[target_taxa_col]] %in% top_taxa, .data[[target_taxa_col]], "Other")) %>%
  mutate(across(c(host_taxa_reg, location_type, vFAM, vGENUS), as.factor))

RANDOM_EFF <- "vFAM"

besthit_dict <- annot %>%
  filter(best_hit != "") %>%
  group_by(best_hit) %>% slice(1) %>% ungroup() %>%
  select(Feature = best_hit, Product = product)

pc_map_linked <- pc_map %>%
  mutate(Phage_Name = str_remove(protein_id, "_CDS_[0-9]+$")) %>%
  inner_join(meta_clean %>% select(Phage_Name, vOTU), by = "Phage_Name")

get_prev_matrix <- function(link_df, otu_df) {
  link_df %>%
    distinct(vOTU, Feature, Phage_Name) %>%
    count(vOTU, Feature, name = "hits") %>%
    inner_join(otu_df %>% select(vOTU, total_genomes), by = "vOTU") %>%
    mutate(prev = hits / total_genomes) %>%
    mutate(prev = if_else(prev > 1, 1, prev)) %>%
    group_by(Feature) %>%
    filter(sum(prev >= 0.4) >= 5) %>%
    ungroup() %>%
    select(vOTU, Feature, prev) %>%
    pivot_wider(names_from = Feature, values_from = prev, values_fill = 0)
}

run_two_step_analysis <- function(feature_mat, otu_df) {
  df <- otu_df %>% inner_join(feature_mat, by = "vOTU")
  features <- setdiff(names(feature_mat), "vOTU")
  fixed_f <- "host_taxa_reg + location_type + log10(total_genomes)"
  
  glm_res <- data.frame()
  for(f in features) {
    form <- paste0("cbind(active_count, inactive_count) ~ `", f, "` + ", fixed_f)
    try({
      m <- glm(as.formula(form), data = df, family = quasibinomial)
      coef <- summary(m)$coefficients
      glm_res <- rbind(glm_res, data.frame(
        Feature = f,
        GLM_Est = coef[2, 1],
        GLM_P = coef[2, 4]
      ))
    }, silent = TRUE)
  }
  
  if(nrow(glm_res) == 0) return(NULL)
  
  glm_res$GLM_FDR <- p.adjust(glm_res$GLM_P, method = "BH")
  candidates <- glm_res %>% filter(GLM_FDR <= 0.05) %>% pull(Feature)
  
  if(length(candidates) == 0) return(NULL)
  
  glmm_res <- data.frame()
  for(f in candidates) {
    form <- paste0("cbind(active_count, inactive_count) ~ `", f, "` + ", fixed_f, " + (1 | ", RANDOM_EFF, ")")
    try({
      m <- glmer(as.formula(form), data = df, family = binomial,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=1e5)))
      coef <- summary(m)$coefficients
      row_idx <- which(grepl(f, rownames(coef), fixed = TRUE))
      if(length(row_idx) == 0) row_idx <- 2
      glmm_res <- rbind(glmm_res, data.frame(
        Feature = f,
        GLMM_Est = coef[row_idx, 1],
        GLMM_P = coef[row_idx, 4],
        Singular = isSingular(m)
      ))
    }, silent = TRUE)
  }
  
  final <- glm_res %>%
    filter(Feature %in% candidates) %>%
    left_join(glmm_res, by = "Feature") %>%
    mutate(GLMM_FDR = p.adjust(GLMM_P, method = "BH")) %>%
    arrange(GLMM_FDR)
  
  return(final)
}

link_annot <- pc_map_linked %>%
  inner_join(annot, by = c("rep_seq_id" = "target")) %>%
  filter(best_hit != "" & !grepl("hypothetical", product, ignore.case = TRUE)) %>%
  select(vOTU, Feature = best_hit, Phage_Name)

mat_annot <- get_prev_matrix(link_annot, otu_stats)
res_annot <- run_two_step_analysis(mat_annot, otu_stats)

res_df <- res_annot %>% left_join(besthit_dict, by = "Feature")
write.csv(res_df, paste0(OUT_PREFIX, "Statistical_Results.csv"), row.names = FALSE)

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

plot_summary_list <- list()
detail_raw_list   <- list()

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

all_details <- bind_rows(detail_raw_list) %>%
  arrange(factor(Feature, levels = target_feats), desc(Group), desc(activity_rate))

write.csv(all_details, paste0(OUT_PREFIX, "Detail_Statistics.csv"), row.names = FALSE)

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

ggsave(paste0(OUT_PREFIX, "Top10_Plot.pdf"), combined_plot, width = 15, height = 7)
ggsave(paste0(OUT_PREFIX, "Top10_Plot.png"), combined_plot, width = 15, height = 7, dpi = 300)
