suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

data_dir     <- "/home/huangyan/experiment/SPA/2-2/used_results"
out_data_tsv <- "phage_correlation_data.tsv"  
out_plot_dir <- "./plot_output/correlation"   
sample_names <- c("Con_AB", "M_AB", "GBIC")
pt_color   <- "#606F7B"   
pt_size    <- 2.5         
pt_alpha   <- 0.7         

line_color <- "black"     
line_size  <- 1.2         
fill_color <- "darkgray"  
fill_alpha <- 0.6         
ref1_val <- 0.47
ref1_col <- "darkgray"
ref1_wid <- 0.8
ref1_typ <- "dashed"   

ref2_val <- 0.70
ref2_col <- "#E57373"
ref2_wid <- 0.8
ref2_typ <- "dashed"
shade_col_gray   <- "#F2F2F2"  # (0, 0.47)
shade_col_blue   <- "#E1F5FE"  # [0.47, 0.70)
shade_col_purple <- "#F3E5F5"  # [0.70, 1.0)
shade_alpha      <- 0.6        
hl_fill   <- "#E41A1C"    
hl_border <- "black"      
hl_size   <- 4.5          
hl_stroke <- 1.2         
target_phages <- data.frame(
  Folder = c("DA914", "DA377", "DA55", "DA14"),
  prophage = c("DA914_16_1_35324", "DA377_11_1_39359", "DA55_29_1_39707", "DA14_9_64261_147917")
)

global_dpi <- 300
plot_w     <- 6
plot_h     <- 6

border_col <- "black"
border_w   <- 1

axis_title_size <- 14
axis_title_col  <- "black"
axis_title_face <- "bold"

axis_text_size  <- 12
axis_text_col   <- "black"
axis_text_face  <- "bold"
axis_text_x_angle <- 0      
axis_text_x_hjust <- 0.5    

axis_tick_col   <- "black"  
axis_tick_w     <- 0.8      

text_size    <- 3.5      
text_v_just  <- -0.3     
title_margin <- 30       

if (!dir.exists(out_plot_dir)) dir.create(out_plot_dir, recursive = TRUE)

stats_files <- list.files(data_dir, pattern = "basic_stats_summary\\.csv$", full.names = TRUE)
if(length(stats_files) == 0) stop("not find basic_stats_summary.csv 文件！")

df_stats <- bind_rows(lapply(stats_files, suppressWarnings(fread)))

host_data <- df_stats %>%
  filter(Region == "host_non_phage") %>%
  select(Folder, Host_Median = Median) %>%
  distinct(Folder, .keep_all = TRUE)

phage_data <- df_stats %>%
  filter(Region == "phage") %>%
  mutate(prophage = gsub("^phage_depth_|\\.txt$", "", File)) %>%
  select(Folder, prophage, Phage_Median = Median) %>%
  distinct(Folder, prophage, .keep_all = TRUE)

df_ratio <- phage_data %>%
  left_join(host_data, by = "Folder") %>%
  mutate(
    Ratio = Phage_Median / Host_Median,
    Log2_Ratio = ifelse(!is.na(Ratio) & Ratio > 0, log2(Ratio), NA)
  ) %>%
  select(Folder, prophage, Log2_Ratio)

list_spa_p <- list()
list_spa_s <- list()

for (sample in sample_names) {
  spa_p_file <- file.path(data_dir, paste0(sample, "_spa_precise.csv"))
  spa_s_file <- file.path(data_dir, paste0(sample, "_spa_sensitive.csv"))
  
  if (file.exists(spa_p_file)) {
    list_spa_p[[sample]] <- suppressWarnings(fread(spa_p_file)) %>%
      mutate(prophage = gsub("phage_depth_|\\.txt", "", Phage_File)) %>%
      select(Folder, prophage, Delta_Mean)
  }
  if (file.exists(spa_s_file)) {
    list_spa_s[[sample]] <- suppressWarnings(fread(spa_s_file)) %>%
      mutate(prophage = gsub("phage_depth_|\\.txt", "", Phage_File)) %>%
      select(Folder, prophage, Effect_Size_Value)
  }
}

df_spa_p_all <- bind_rows(list_spa_p)
df_spa_s_all <- bind_rows(list_spa_s)

raw_data <- df_ratio %>%
  full_join(df_spa_p_all, by = c("Folder", "prophage")) %>%
  full_join(df_spa_s_all, by = c("Folder", "prophage")) %>%
  mutate(
    Delta_Mean = as.numeric(Delta_Mean),
    Effect_Size_Value = as.numeric(Effect_Size_Value)
  )

write.table(raw_data, file.path(out_plot_dir, out_data_tsv), sep = "\t", row.names = FALSE, quote = FALSE)


plot_cor_final <- function(data, x_col, y_col, x_lab, y_lab, title_text, filename_out, 
                           add_lines=FALSE, add_shade=FALSE, hl_targets=NULL) {

  df <- data.frame(
    Folder = data$Folder,
    prophage = data$prophage,
    x = data[[x_col]], 
    y = data[[y_col]]
  )
  df <- df[!is.na(df$x) & !is.na(df$y), ]

  df_hl <- data.frame()
  if (!is.null(hl_targets)) {
    df_hl <- df %>% inner_join(hl_targets, by = c("Folder", "prophage"))
  }

  if (nrow(df) > 2) {
    p_cor <- cor.test(df$x, df$y, method = "pearson")
    s_cor <- suppressWarnings(cor.test(df$x, df$y, method = "spearman")) 

    fmt_p <- function(p) if(p < 0.001) "< 0.001" else sprintf("= %.3f", p)
    label_text <- paste0(
      "Pearson: R = ", round(p_cor$estimate, 3), ", p ", fmt_p(p_cor$p.value), "\n",
      "Spearman: R = ", round(s_cor$estimate, 3), ", p ", fmt_p(s_cor$p.value)
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      {if(add_shade) list(
        annotate("rect", xmin = -Inf, xmax = ref1_val, ymin = -Inf, ymax = ref1_val, fill = shade_col_gray, alpha = shade_alpha),
        annotate("rect", xmin = ref1_val, xmax = ref2_val, ymin = ref1_val, ymax = ref2_val, fill = shade_col_blue, alpha = shade_alpha),
        annotate("rect", xmin = ref2_val, xmax = Inf, ymin = ref2_val, ymax = Inf, fill = shade_col_purple, alpha = shade_alpha)
      )} +
      geom_point(color = pt_color, size = pt_size, alpha = pt_alpha) +
      geom_smooth(method = "lm", color = line_color, fill = fill_color, alpha = fill_alpha, linewidth = line_size) +
      {if(add_lines) list(
        geom_vline(xintercept = ref1_val, color = ref1_col, linewidth = ref1_wid, linetype = ref1_typ),
        geom_hline(yintercept = ref1_val, color = ref1_col, linewidth = ref1_wid, linetype = ref1_typ),
        geom_vline(xintercept = ref2_val, color = ref2_col, linewidth = ref2_wid, linetype = ref2_typ),
        geom_hline(yintercept = ref2_val, color = ref2_col, linewidth = ref2_wid, linetype = ref2_typ)
      )} +
      {if(nrow(df_hl) > 0) list(
        geom_point(data = df_hl, aes(x = x, y = y), 
                   fill = hl_fill, color = hl_border, size = hl_size, shape = 21, stroke = hl_stroke)
      )} +

      labs(title = title_text, x = x_lab, y = y_lab) +
      theme_classic(base_family = "Arial") +

      theme(
        panel.border = element_rect(colour = border_col, fill = NA, linewidth = border_w),
        panel.grid = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_line(color = axis_tick_col, linewidth = axis_tick_w),
        axis.title.x = element_text(size = axis_title_size, color = axis_title_col, face = axis_title_face),
        axis.title.y = element_text(size = axis_title_size, color = axis_title_col, face = axis_title_face),
        axis.text.x = element_text(size = axis_text_size, color = axis_text_col, face = axis_text_face, angle = axis_text_x_angle, hjust = axis_text_x_hjust),
        axis.text.y = element_text(size = axis_text_size, color = axis_text_col, face = axis_text_face),
        plot.margin = margin(t = 50, r = 20, b = 20, l = 20, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = title_margin))
      ) +
      coord_cartesian(clip = "off")

    ggsave(file.path(out_plot_dir, filename_out), p, width = plot_w, height = plot_h, dpi = global_dpi, device = cairo_pdf)
    
 
    if (!is.null(hl_targets)) {
      message(sprintf("       -> highlight %d phage！", nrow(df_hl)))
    }

  } else {
    message("     error")
  }
}

# Delta_Mean vs Effect_Size_Value
plot_cor_final(
  data = raw_data, 
  x_col = "Delta_Mean", y_col = "Effect_Size_Value",
  x_lab = "SPA Precise: Delta Mean", y_lab = "SPA Sensitive: Effect Size",
  title_text = "Correlation: Delta Mean vs Effect Size",
  filename_out = "04_cor_delta_effect_final.pdf",
  add_lines = TRUE,
  add_shade = TRUE,               
  hl_targets = target_phages      
)


