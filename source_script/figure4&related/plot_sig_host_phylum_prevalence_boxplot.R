library(tidyverse)
library(ggpubr)

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

generate_prevalence_boxplot <- function(input_file, level_name, out_file) {
  if (!file.exists(input_file)) {
    message("Warning: Input file not found - ", input_file)
    return(NULL)
  }
  
  message("\n>>> Processing data for: ", level_name)
  message("Reading data from: ", input_file)

  stats_df <- read_tsv(input_file, show_col_types = FALSE)

  sig_df <- stats_df %>%
    filter(Significance_Type == "Significant") %>%
    group_by(Phylum) %>%
    mutate(Phylum_Count = n()) %>%
    ungroup() %>%
    mutate(Phylum_Label = paste0(Phylum, "\n(n=", Phylum_Count, ")")) %>%
    mutate(Phylum_Label = fct_reorder(Phylum_Label, Prevalence_Pct, .fun = median, .desc = TRUE))
  
  message("Total Exact-Significant ", level_name, "s filtered: ", nrow(sig_df))
  
  if (nrow(sig_df) == 0) {
    message("No significant entities found. Skipping plot generation for ", level_name)
    return(NULL)
  }
  
  actual_phyla <- unique(sig_df$Phylum)
  current_palette <- c()
  for(phylum in actual_phyla) {
    if(phylum %in% names(PHYLUM_COLORS)) {
      current_palette[phylum] <- PHYLUM_COLORS[phylum]
    } else {
      current_palette[phylum] <- "#bdbdbd"
    }
  }
  
  p_box <- ggplot(sig_df, aes(x = Phylum_Label, y = Prevalence_Pct, fill = Phylum)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6, color = "black") +
    geom_jitter(aes(size = Genome_Count),
                width = 0.2, alpha = 0.9, shape = 21, color = "black", stroke = 0.6) +
    scale_fill_manual(values = current_palette) +
    scale_size_continuous(range = c(2, 8), name = "Genome Count") +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, by = 20)) +
    labs(
      title = paste0("Prevalence Distribution of Significantly Active ", level_name, "s"),
      subtitle = paste0("Total Exact-Significant ", level_name, "s: ", nrow(sig_df)),
      x = "Host Phylum",
      y = "Active Prevalence (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      plot.title = element_text(size = 16, face = "bold", color = "black"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.minor = element_blank()
    ) +
    guides(fill = "none", size = guide_legend(override.aes = list(fill = "gray50"))) +
    stat_compare_means(method = "kruskal.test", label.y = 105, size = 4.5, color = "red3", fontface = "italic")
  
  ggsave(out_file, p_box, width = 9, height = 7, bg = "white")
  message("Boxplot successfully saved to: ", out_file)
}

generate_prevalence_boxplot(
  input_file = "stats_baseline_vOTU_EXACT.tsv",
  level_name = "vOTU",
  out_file   = "Prevalence_Boxplot_Significant_vOTU.pdf"
)

generate_prevalence_boxplot(
  input_file = "stats_baseline_vFAM_EXACT.tsv",
  level_name = "vFAM",
  out_file   = "Prevalence_Boxplot_Significant_vFAM.pdf"
)

message("\n>>> All tasks completed successfully.")
