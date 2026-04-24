
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggsci)
  library(scales)
})

file_path    <- "../prophage_stats_final.tsv"
out_dir      <- "./plot_output"  

plot_6panel_w      <- 14       
plot_6panel_h      <- 15       

plot_method_w      <- 12
plot_method_h      <- 7

global_dpi         <- 600     
global_base_family <- "Arial"  

legend_position    <- "top"   
legend_text_size   <- 14      
legend_text_color  <- "black"  

title_size         <- 12      
title_face         <- "bold"   

color_axis_title   <- "black"
size_axis_title    <- 18
face_axis_title    <- "bold"

color_axis_text    <- "black"
size_axis_text     <- 16
face_axis_text     <- "plain"

axis_text_angle_x  <- 45       
panel_border_color <- "black"  
panel_border_width <- 1       
color_lysogeny     <- "#374E55FF" 
color_induction    <- "#DF8F44FF" 
color_poly         <- "#00A1D5FF" 
metric_colors      <- setNames(c(color_lysogeny, color_induction, color_poly), 
                               c("Lysogeny_Rate", "Induction_Rate", "Poly_Lysogeny_Rate"))

bar_width          <- 0.7
bar_border_color   <- "black"
bar_border_width   <- 0.5
bar_label_size     <- 4      
bar_label_face     <- "bold"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df <- read_tsv(file_path, show_col_types = FALSE) %>%
  mutate(
    Prophage_Count = replace_na(as.numeric(Prophage_Count), 0),
    Active_Prophage_Count = replace_na(as.numeric(Active_Prophage_Count), 0)
  )

df_total <- df %>% mutate(source = "Total")
df_combined <- bind_rows(df, df_total)
source_levels <- c(sort(unique(df$source)), "Total")
df_combined$source <- factor(df_combined$source, levels = source_levels)

metric_levels <- c("Lysogeny_Rate", "Induction_Rate", "Poly_Lysogeny_Rate")

shared_theme <- theme_bw(base_size = 12, base_family = global_base_family) +
  theme(
    plot.title   = element_text(hjust = 0, face = title_face, size = title_size, color = "black"),
    
    axis.title.x = element_text(size = size_axis_title, color = color_axis_title, face = face_axis_title),
    axis.title.y = element_text(size = size_axis_title, color = color_axis_title, face = face_axis_title),
    
    axis.text.x  = element_text(size = size_axis_text, color = color_axis_text, face = face_axis_text, angle = axis_text_angle_x, hjust = ifelse(axis_text_angle_x > 0, 1, 0.5)),
    axis.text.y  = element_text(size = size_axis_text, color = color_axis_text, face = face_axis_text),
    
    panel.grid   = element_blank(),
    panel.border = element_rect(colour = panel_border_color, fill = NA, linewidth = panel_border_width),
    
    legend.position   = legend_position,
    legend.title      = element_blank(),
    legend.text       = element_text(size = legend_text_size, color = legend_text_color),
    legend.background = element_blank()
  )

calculate_stats <- function(data, tax_rank = NULL, group_col = "source") {
  
  if (is.null(tax_rank) || tax_rank == "Unweighted") {
    stats <- data %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        Lysogeny_Rate = mean(Prophage_Count > 0),
        n_lysogens = sum(Prophage_Count > 0),
        n_active = sum(Active_Prophage_Count > 0),
        n_poly_events = sum(Prophage_Count >= 2),
        Induction_Rate = if_else(n_lysogens > 0, n_active / n_lysogens, 0),
        Poly_Lysogeny_Rate = if_else(n_lysogens > 0, n_poly_events / n_lysogens, 0),
        .groups = "drop"
      )
  } else {
    taxa_stats <- data %>%
      group_by(.data[[group_col]], .data[[tax_rank]]) %>%
      summarise(
        t_lysogeny = mean(Prophage_Count > 0),
        t_n_lysogens = sum(Prophage_Count > 0),
        t_n_active = sum(Active_Prophage_Count > 0),
        t_n_poly = sum(Prophage_Count >= 2),
        t_induction = if_else(t_n_lysogens > 0, t_n_active / t_n_lysogens, NA_real_),
        t_poly = if_else(t_n_lysogens > 0, t_n_poly / t_n_lysogens, NA_real_),
        .groups = "drop"
      )
    
    stats <- taxa_stats %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        Lysogeny_Rate = mean(t_lysogeny, na.rm = TRUE),
        Induction_Rate = mean(t_induction, na.rm = TRUE),
        Poly_Lysogeny_Rate = mean(t_poly, na.rm = TRUE),
        .groups = "drop"
      )
  }
  return(stats)
}


plot_source_comparison <- function(data, tax_rank, title_suffix) {
  
  stats_df <- calculate_stats(data, tax_rank = tax_rank, group_col = "source")
  
  plot_data <- stats_df %>%
    pivot_longer(cols = ends_with("_Rate"), names_to = "Metric", values_to = "Rate") %>%
    mutate(Metric = factor(Metric, levels = metric_levels)) 
  
  p <- ggplot(plot_data, aes(x = source, y = Rate, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.8), width = bar_width, 
             color = bar_border_color, linewidth = bar_border_width) +
    geom_text(aes(label = sprintf("%.2f", Rate)),
              position = position_dodge(width = 0.8), vjust = -0.5, 
              size = bar_label_size, fontface = bar_label_face) +
    scale_fill_manual(values = metric_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 1)) +
    labs(title = paste(title_suffix), y = "Rate", x = NULL) +
    shared_theme
  
  return(p)
}

p1 <- plot_source_comparison(df_combined, NULL, "1. Unweighted (Original)")
p2 <- plot_source_comparison(df_combined, "Phylum", "2. Weighted by Phylum")
p3 <- plot_source_comparison(df_combined, "Class", "3. Weighted by Class")
p4 <- plot_source_comparison(df_combined, "Order", "4. Weighted by Order")
p5 <- plot_source_comparison(df_combined, "Family", "5. Weighted by Family")
p6 <- plot_source_comparison(df_combined, "Genus", "6. Weighted by Genus")

final_design_6panel <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_layout(guides = "collect") & theme(legend.position = legend_position)

out_file_1 <- file.path(out_dir, "Figure_Weighted_Source_Comparison.pdf")

ggsave(out_file_1, final_design_6panel, width = plot_6panel_w, height = plot_6panel_h, dpi = global_dpi, device = cairo_pdf)

rank_levels <- c("Unweighted", "Phylum", "Class", "Order", "Family", "Genus", "Species")
list_results <- list()

for (r in rank_levels) {
  temp_df <- df %>% mutate(dummy_group = "All_Data")
  res <- calculate_stats(temp_df, tax_rank = if(r == "Unweighted") NULL else r, group_col = "dummy_group")
  res$Method <- r 
  list_results[[r]] <- res
}

df_methods <- bind_rows(list_results) %>%
  select(-dummy_group) %>%
  pivot_longer(cols = ends_with("_Rate"), names_to = "Metric", values_to = "Rate") %>%
  mutate(
    Metric = factor(Metric, levels = metric_levels), 
    Method = factor(Method, levels = rank_levels)    
  )

p_method <- ggplot(df_methods, aes(x = Method, y = Rate, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.8), width = bar_width, 
           color = bar_border_color, linewidth = bar_border_width) +
  geom_text(aes(label = sprintf("%.2f", Rate)),
            position = position_dodge(width = 0.8), vjust = -0.5, 
            size = bar_label_size + 0.5, fontface = bar_label_face) +
  scale_fill_manual(values = metric_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 1)) +
  labs(
    title = "Impact of Taxonomic Weighting on Prophage Rates (Total Dataset)",
    subtitle = "Definition: Induction Rate = Active Lysogens / All Lysogens",
    x = "Weighting Granularity",
    y = "Calculated Rate"
  ) +
  shared_theme

out_file_2 <- file.path(out_dir, "Figure_Total_Method_Comparison.pdf")

ggsave(out_file_2, p_method, width = plot_method_w, height = plot_method_h, dpi = global_dpi, device = cairo_pdf)

