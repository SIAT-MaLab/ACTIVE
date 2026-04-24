
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggdist)
  library(ggpubr)
  library(dplyr)
})


input_file <- "phage_full_stats_separated.tsv"

raw <- read.delim(input_file, stringsAsFactors=FALSE, check.names = FALSE)

raw$Activity_Score <- as.numeric(as.character(raw$Activity_Score))
raw$Log2_Ratio <- as.numeric(as.character(raw$Log2_Ratio))

raw$Activity_Score[is.na(raw$Activity_Score)] <- 0

clean_df <- raw[!is.na(raw$Log2_Ratio) & is.finite(raw$Log2_Ratio), ]

draw_raincloud_horizontal <- function(data_vec, title_text) {

  df_plot <- data.frame(value = data_vec)

  stats <- quantile(df_plot$value, probs = c(0.25, 0.5, 0.75))
  q1 <- stats[1]
  median_val <- stats[2]
  q3 <- stats[3]

  cat(sprintf("   Median: %.4f, Q1: %.4f, Q3: %.4f\n", median_val, q1, q3))
  
  p <- ggplot(df_plot, aes(x = value, y = 1)) +
    stat_halfeye(
      adjust = 0.6, 
      justification = -0.2, 
      .width = 0, 
      point_colour = NA,
      fill = "#E0E0E0",        
      slab_color = "black",    
      slab_linewidth = 0.5, 
      alpha = 1.0
    ) +
    geom_jitter(
      aes(y = 1), 
      height = 0.08, 
      size = 1.2, 
      alpha = 0.6, 
      color = "#4D4D4D"        
    ) +

    geom_vline(xintercept = median_val, color = "#d95f02", size = 1, linetype = "solid") +
    geom_vline(xintercept = c(q1, q3), color = "#333333", size = 0.6, linetype = "dashed") +

    labs(x = "Log2 Ratio", y = "", title = title_text) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.line.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(size = 16),
      axis.text.x  = element_text(size = 14),
    ) +
    scale_x_continuous(n.breaks = 8)
  
  return(p)
}

vec_high <- clean_df$Log2_Ratio[clean_df$Activity_Score >= 0.7]
p_high <- draw_raincloud_horizontal(vec_high, "High Activity (Score ≥ 0.7)")

vec_all <- clean_df$Log2_Ratio
p_all <- draw_raincloud_horizontal(vec_all, "All Data (Full Set)")

final_plot <- ggarrange(p_high, p_all, 
                        ncol = 1, 
                        nrow = 2,
                        heights = c(1, 1)) 

output_filename <- "raincloud_horizontal_comparison.pdf"

ggsave(output_filename, final_plot, width = 8, height = 5)

