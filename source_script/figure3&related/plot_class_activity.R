
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

input_file         <- "final_merged_phage_metadata.tsv"
out_file_all       <- "class_activity.pdf"
out_file_hq        <- "hq_class_activity.pdf"

out_tsv_all        <- "class_activity_stats.tsv"
out_tsv_hq         <- "hq_class_activity_stats.tsv"

threshold_activity     <- 0.7
threshold_completeness <- 90

plot_w             <- 11
plot_h             <- 6

color_total        <- "#3b262b"
color_target       <- "#A41A00"
color_normal       <- "#AEAFB0"
color_other        <- "#F2F3F4"

width_bar          <- 0.7
color_bar_border   <- "black"
width_bar_border   <- 0.3

size_title         <- 14
size_axis_title    <- 18
size_axis_text     <- 16
size_bar_label     <- 4

color_axis_title    <- "black"
color_axis_text_nor <- "grey40"
color_axis_text_hl  <- "black"
color_bar_label     <- "black"


get_summary_stats <- function(data) {
  summary_grouped <- data %>%
    group_by(ictv_Class) %>%
    summarise(
      Counts = n(),
      Active_Count = sum(is_active, na.rm = TRUE),
      Active_Rate = Active_Count / Counts,
      .groups = "drop"
    )

  summary_total <- data %>%
    summarise(
      ictv_Class = "Total",
      Counts = n(),
      Active_Count = sum(is_active, na.rm = TRUE),
      Active_Rate = Active_Count / Counts
    )


  summary_df <- bind_rows(summary_grouped, summary_total) %>%
    mutate(
      Type = case_when(
        ictv_Class == "Total" ~ "Total",
        ictv_Class == "Caudoviricetes" ~ "Target",
        ictv_Class == "Other" ~ "Other",
        TRUE ~ "Normal"
      )
    ) %>%
    arrange(Active_Rate)
  
  return(summary_df)
}

generate_activity_plot <- function(summary_df) {

  summary_df <- summary_df %>%
    mutate(ictv_Class = factor(ictv_Class, levels = c(setdiff(ictv_Class, "Total"), "Total")))

  y_axis_levels <- levels(summary_df$ictv_Class)
  color_mapping <- c("Target" = color_target, "Normal" = color_normal, "Other" = color_other, "Total" = color_total)

  pA <- ggplot(summary_df, aes(x = Active_Rate, y = ictv_Class, fill = Type)) +
    geom_col(width = width_bar, show.legend = FALSE, color = color_bar_border, linewidth = width_bar_border) +
    geom_text(aes(label = percent(Active_Rate, accuracy = 0.1)),
              hjust = -0.1, size = size_bar_label, color = color_bar_label) +
    scale_fill_manual(values = color_mapping) +
    scale_x_continuous(labels = percent, expand = expansion(mult = c(0, 0.2))) +
    labs(title = "A. Class Activity", x = "Active Rate", y = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title   = element_text(hjust = 0, size = size_title, face = "bold"),
      axis.title.x = element_text(size = size_axis_title, color = color_axis_title),
      axis.text.x  = element_text(size = size_axis_text, color = color_axis_text_nor),
      axis.text.y  = element_text(
        size = size_axis_text,
        face = ifelse(y_axis_levels %in% c("Total", "Caudoviricetes"), "bold", "plain"),
        color = ifelse(y_axis_levels %in% c("Total", "Caudoviricetes"), color_axis_text_hl, color_axis_text_nor)
      ),
      axis.ticks   = element_blank()
    )

  pB <- ggplot(summary_df, aes(x = Counts, y = ictv_Class, fill = Type)) +
    geom_col(width = width_bar, show.legend = FALSE, color = color_bar_border, linewidth = width_bar_border) +
    geom_text(aes(label = comma(Counts)),
              hjust = -0.1, size = size_bar_label, color = color_bar_label) +
    scale_fill_manual(values = color_mapping) +
    scale_x_log10(labels = comma, expand = expansion(mult = c(0, 0.2))) +
    labs(title = "B. Genome Counts (Log Scale)", x = "Counts", y = NULL) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title   = element_text(hjust = 0, size = size_title, face = "bold"),
      axis.title.x = element_text(size = size_axis_title, color = color_axis_title),
      axis.text.x  = element_text(size = size_axis_text, color = color_axis_text_nor),
      axis.text.y  = element_blank(),
      axis.ticks   = element_blank()
    )

  return(pA + pB + plot_layout(widths = c(1, 1)))
}


df_clean <- read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(Activity_Score = as.numeric(Activity_Score)) %>%
  mutate(
    is_active = ifelse(!is.na(Activity_Score) & Activity_Score >= threshold_activity, 1, 0),
    ictv_Class = ifelse(is.na(ictv_Class) | ictv_Class == "", "Other", ictv_Class),
    ictv_Class = str_remove(ictv_Class, "^[a-z]__")
  )

stats_all <- get_summary_stats(df_clean)
write_tsv(stats_all, out_tsv_all) 
plot_all <- generate_activity_plot(stats_all)
ggsave(out_file_all, plot = plot_all, width = plot_w, height = plot_h, device = "pdf")

df_hq <- df_clean %>%
  filter(!is.na(completeness)) %>%
  filter(as.numeric(completeness) >= threshold_completeness)

if (nrow(df_hq) > 0) {
  stats_hq <- get_summary_stats(df_hq)
  write_tsv(stats_hq, out_tsv_hq) 
  plot_hq <- generate_activity_plot(stats_hq)
  ggsave(out_file_hq, plot = plot_hq, width = plot_w, height = plot_h, device = "pdf")
} else {
  message("withou hq")
}

