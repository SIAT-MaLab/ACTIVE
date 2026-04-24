suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

data_dir <- "/home/huangyan/experiment/SPA/2-2/used_results"
out_dir  <- "./plot_output/activity_composition"
sample_names <- c("Con_AB", "M_AB", "GBIC")

spa_threshold_high <- 0.7
spa_threshold_low  <- 0.47

label_sig_active <- "Significant Active"
label_amb_active <- "Ambiguous Active"
label_dormant    <- "Dormant"
label_unassess   <- "Unassessable"

status_colors <- setNames(
  c("#E41A1C", "#FF7F00", "#377EB8", "grey70"),
  c(label_sig_active, label_amb_active, label_dormant, label_unassess)
)

global_dpi <- 300
plot_w     <- 6
plot_h     <- 5
bar_border_col <- "black"   
bar_border_w   <- 0.3        
bar_alpha      <- 0.9        
axis_title_size <- 14
axis_title_col  <- "black"
axis_title_face <- "bold"    # "plain", "bold", "italic"
axis_text_size  <- 12
axis_text_col   <- "black"
axis_text_face  <- "bold"
axis_text_x_angle <- 45      
axis_text_x_hjust <- 1      
axis_line_col   <- "black"   
axis_line_w     <- 0.8       
axis_tick_col   <- "black"   
axis_tick_w     <- 0.8       
legend_pos        <- "right" 
legend_text_size  <- 11
legend_text_col   <- "black"
legend_text_face  <- "plain"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


for (sample in sample_names) {
  prop_file  <- file.path(data_dir, paste0(sample, "_propagate.tsv"))
  spa_p_file <- file.path(data_dir, paste0(sample, "_spa_precise.csv"))
  spa_s_file <- file.path(data_dir, paste0(sample, "_spa_sensitive.csv"))
  
  df_prop <- suppressWarnings(fread(prop_file)) %>%
    filter(rowSums(is.na(select(., -prophage))) < (ncol(.) - 1)) %>%
    mutate(
      Status = case_when(
        host_len == 0 ~ label_unassess,
        active == "active" ~ label_sig_active,
        active == "dormant" ~ label_dormant,
        TRUE ~ label_unassess
      ),
      Method = "PropagAtE"
    ) %>% select(prophage, Method, Status)
  

  df_spa_p <- suppressWarnings(fread(spa_p_file)) %>%
    mutate(
      prophage = gsub("phage_depth_|\\.txt", "", Phage_File),
      Value = as.numeric(Delta_Mean),
      Status = case_when(
        is.na(Value) ~ label_dormant,
        Value >= spa_threshold_high ~ label_sig_active,
        Value >= spa_threshold_low & Value < spa_threshold_high ~ label_amb_active,
        TRUE ~ label_dormant
      ),
      Method = "SPA Precise"
    ) %>% select(prophage, Method, Status)
  

  df_spa_s <- suppressWarnings(fread(spa_s_file)) %>%
    mutate(
      prophage = gsub("phage_depth_|\\.txt", "", Phage_File),
      Value = as.numeric(Effect_Size_Value),
      Status = case_when(
        is.na(Value) ~ label_dormant,
        Value >= spa_threshold_high ~ label_sig_active,
        Value >= spa_threshold_low & Value < spa_threshold_high ~ label_amb_active,
        TRUE ~ label_dormant
      ),
      Method = "SPA Sensitive"
    ) %>% select(prophage, Method, Status)
  

  df_all <- bind_rows(df_prop, df_spa_p, df_spa_s)
  
  df_summary <- df_all %>%
    group_by(Method, Status) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Method = factor(Method, levels = c("PropagAtE", "SPA Precise", "SPA Sensitive")),
      Status = factor(Status, levels = c(label_sig_active, label_amb_active, label_dormant, label_unassess))
    )
  
  df_total <- df_summary %>%
    group_by(Method) %>%
    summarise(TotalCount = sum(Count))
  

  p <- ggplot(df_summary, aes(x = Method, y = Count, fill = Status)) +
    geom_col(color = bar_border_col, linewidth = bar_border_w, alpha = bar_alpha, width = 0.6) +
    scale_fill_manual(values = status_colors) +
    theme_classic(base_family = "Arial") +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(size = axis_title_size, color = axis_title_col, face = axis_title_face),
      axis.title.y = element_text(size = axis_title_size, color = axis_title_col, face = axis_title_face),
      axis.text.x  = element_text(size = axis_text_size, color = axis_text_col, face = axis_text_face, 
                                  angle = axis_text_x_angle, hjust = axis_text_x_hjust),
      axis.text.y  = element_text(size = axis_text_size, color = axis_text_col, face = axis_text_face),
      axis.line    = element_line(color = axis_line_col, linewidth = axis_line_w),
      axis.ticks   = element_line(color = axis_tick_col, linewidth = axis_tick_w),
      legend.position = legend_pos,
      legend.title    = element_blank(),
      legend.text     = element_text(size = legend_text_size, color = legend_text_col, face = legend_text_face)
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = paste0("Activity Classification: ", sample), x = "", y = "Number of Prophages")
  out_name <- paste0(sample, "_composition.pdf")
  ggsave(file.path(out_dir, out_name), p, width = plot_w, height = plot_h, dpi = global_dpi, device = cairo_pdf)

}


