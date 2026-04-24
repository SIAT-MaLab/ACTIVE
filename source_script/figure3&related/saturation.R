#!/usr/bin/env Rscript


rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales) 
})

output_dir <- "saturation"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

custom_phylum_colors <- c(
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


main_file <- "merged_phage_stats_taxonomy.tsv"

input_df <- read.delim(main_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "")
input_df$Activity_Score <- as.numeric(as.character(input_df$Activity_Score))

paths_votu <- c("allpp_taxonomy/active_tax/stats_fisher_vOTU.tsv", "stats_fisher_vOTU.tsv")
stats_file_votu <- NULL
for (p in paths_votu) { if (file.exists(p)) { stats_file_votu <- p; break } }
stats_df_votu <- read.delim(stats_file_votu, header = TRUE, sep = "\t", check.names = FALSE, quote = "")

# vFAM
paths_vfam <- c("allpp_taxonomy/active_tax/stats_fisher_vFAM.tsv", "stats_fisher_vFAM.tsv", 
                "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_taxonomy/active_tax/stats_fisher_vFAM.tsv")
stats_file_vfam <- NULL
for (p in paths_vfam) { if (file.exists(p)) { stats_file_vfam <- p; break } }
stats_df_vfam <- read.delim(stats_file_vfam, header = TRUE, sep = "\t", check.names = FALSE, quote = "")
list_active_votu <- input_df %>% filter(Activity_Score >= 0.7) %>% pull(vOTU) %>% unique()
list_active_vfam <- input_df %>% filter(Activity_Score >= 0.7) %>% pull(vFAM) %>% unique()
list_sig_votu <- stats_df_votu %>% filter(Significance == "Significant") %>% pull(vOTU) %>% unique()
list_sig_vfam <- stats_df_vfam %>% filter(Significance == "Significant") %>% pull(vFAM) %>% unique()
processed_df <- input_df %>%
  select(Genome, Phylum, vOTU, vFAM) %>% 
  mutate(
    Is_Active_vOTU = vOTU %in% list_active_votu,
    Is_Significant_vOTU = vOTU %in% list_sig_votu,
    Is_Active_vFAM = vFAM %in% list_active_vfam,
    Is_Significant_vFAM = vFAM %in% list_sig_vfam
  ) %>%
  filter(!is.na(Phylum) & Phylum != "")

===============================

calculate_accumulation <- function(data, entity_col, metric_col = NULL, steps = 30, iterations = 10) {

  all_genomes <- unique(data$Genome)
  n_total <- length(all_genomes)

  if (n_total < steps) {
    sample_sizes <- 1:n_total
  } else {
    sample_sizes <- unique(round(seq(1, n_total, length.out = steps)))
    if (tail(sample_sizes, 1) != n_total) sample_sizes <- c(sample_sizes, n_total)
  }

  target_data <- data %>% filter(!is.na(.data[[entity_col]]) & .data[[entity_col]] != "")
  
  if (!is.null(metric_col)) {
    target_data <- target_data %>% filter(.data[[metric_col]] == TRUE)
  }

  results <- data.frame(
    Sample_Depth = integer(),
    Mean_Count = numeric(),
    SD_Count = numeric()
  )

  for (size in sample_sizes) {
    counts <- numeric(iterations)
    for (i in 1:iterations) {
      picked_genomes <- sample(all_genomes, size)
      n_unique <- length(unique(target_data[[entity_col]][target_data$Genome %in% picked_genomes]))
      counts[i] <- n_unique
    }

    results <- rbind(results, data.frame(
      Sample_Depth = size,
      Mean_Count = mean(counts),
      SD_Count = sd(counts)
    ))
  }
  return(results)
}


phylum_stats <- processed_df %>%
  group_by(Phylum) %>%
  summarise(n_genomes = n_distinct(Genome)) %>%
  arrange(desc(n_genomes)) %>%
  filter(n_genomes >= 5) 

target_phyla <- phylum_stats$Phylum
res_list <- list(
  vOTU_total = list(), vOTU_active = list(), vOTU_sig = list(),
  vFAM_total = list(), vFAM_active = list(), vFAM_sig = list()
)

pb <- txtProgressBar(min = 0, max = length(target_phyla), style = 3)

for (i in seq_along(target_phyla)) {
  p_name <- target_phyla[i]
  sub_df <- processed_df %>% filter(Phylum == p_name)


  tmp <- calculate_accumulation(sub_df, "vOTU", NULL); tmp$Phylum <- p_name; res_list$vOTU_total[[i]] <- tmp
  tmp <- calculate_accumulation(sub_df, "vOTU", "Is_Active_vOTU"); tmp$Phylum <- p_name; res_list$vOTU_active[[i]] <- tmp
  tmp <- calculate_accumulation(sub_df, "vOTU", "Is_Significant_vOTU"); tmp$Phylum <- p_name; res_list$vOTU_sig[[i]] <- tmp


  tmp <- calculate_accumulation(sub_df, "vFAM", NULL); tmp$Phylum <- p_name; res_list$vFAM_total[[i]] <- tmp
  tmp <- calculate_accumulation(sub_df, "vFAM", "Is_Active_vFAM"); tmp$Phylum <- p_name; res_list$vFAM_active[[i]] <- tmp
  tmp <- calculate_accumulation(sub_df, "vFAM", "Is_Significant_vFAM"); tmp$Phylum <- p_name; res_list$vFAM_sig[[i]] <- tmp

  setTxtProgressBar(pb, i)
}
close(pb)

df_votu_total  <- do.call(rbind, res_list$vOTU_total)
df_votu_active <- do.call(rbind, res_list$vOTU_active)
df_votu_sig    <- do.call(rbind, res_list$vOTU_sig)
df_vfam_total  <- do.call(rbind, res_list$vFAM_total)
df_vfam_active <- do.call(rbind, res_list$vFAM_active)
df_vfam_sig    <- do.call(rbind, res_list$vFAM_sig)

normalize_curve_data <- function(df) {
  df %>%
    group_by(Phylum) %>%
    mutate(
      Norm_Sample_Depth = Sample_Depth / max(Sample_Depth, na.rm = TRUE),
      Norm_Mean_Count = Mean_Count / max(Mean_Count, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(Norm_Mean_Count = ifelse(is.nan(Norm_Mean_Count), 0, Norm_Mean_Count))
}

df_votu_total  <- normalize_curve_data(df_votu_total)
df_votu_active <- normalize_curve_data(df_votu_active)
df_votu_sig    <- normalize_curve_data(df_votu_sig)

df_vfam_total  <- normalize_curve_data(df_vfam_total)
df_vfam_active <- normalize_curve_data(df_vfam_active)
df_vfam_sig    <- normalize_curve_data(df_vfam_sig)


plot_saturation <- function(data, title_str, y_label_str, file_name) {

  p <- ggplot(data, aes(x = Norm_Sample_Depth, y = Norm_Mean_Count, color = Phylum, group = Phylum)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
    scale_color_manual(values = custom_phylum_colors) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) + 
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = title_str,
      subtitle = "Normalized Saturation curve (Fitted Smooth Lines)",
      x = "Normalized Host Genomes Sampled (%)",
      y = y_label_str
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),                                      
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
      axis.line = element_line(color = "black")                          
    )

  ggsave(paste0(output_dir, "/", file_name), p, width = 10, height = 7)
  
}


plot_saturation(df_votu_total,  "Saturation Curve: Total vOTUs by Phylum",              "Normalized Cumulative vOTUs (%)", "Curve_1_Total_vOTU.pdf")
plot_saturation(df_votu_active, "Saturation Curve: Active vOTUs by Phylum",             "Normalized Cumulative vOTUs (%)", "Curve_2_Active_vOTU.pdf")
plot_saturation(df_votu_sig,    "Saturation Curve: Significant Active vOTUs by Phylum", "Normalized Cumulative vOTUs (%)", "Curve_3_Significant_vOTU.pdf")

plot_saturation(df_vfam_total,  "Saturation Curve: Total vFAMs by Phylum",              "Normalized Cumulative vFAMs (%)", "Curve_4_Total_vFAM.pdf")
plot_saturation(df_vfam_active, "Saturation Curve: Active vFAMs by Phylum",             "Normalized Cumulative vFAMs (%)", "Curve_5_Active_vFAM.pdf")
plot_saturation(df_vfam_sig,    "Saturation Curve: Significant Active vFAMs by Phylum", "Normalized Cumulative vFAMs (%)", "Curve_6_Significant_vFAM.pdf")


