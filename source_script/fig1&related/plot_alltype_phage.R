suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(patchwork) 
  library(parallel)
})

csv_dir <- "/home/huangyan/experiment/SPA/2-2/used_results"
ab_depth_base   <- "/home/huangyan/experiment/SPA/2-2/AB/CM_results/spa_results_20260305/mapping_results"
gbic_depth_base <- "/home/huangyan/experiment/SPA/2-2/GBIC/spa_results_20260303/mapping_results"
out_dir_combined <- "./plot_output/classified_dashboards"
num_threads      <- 20

plot_combined_w    <- 14       
plot_combined_h    <- 6        
global_dpi         <- 300      
global_base_family <- "Arial"   

legend_position    <- "top"    
legend_text_size   <- 16       
legend_text_color  <- "black"  
legend_text_face   <- "plain"  

color_axis_title   <- "black"
size_axis_title    <- 18
face_axis_title    <- "bold"

color_axis_text    <- "black"
size_axis_text     <- 18
face_axis_text     <- "plain"

axis_line_color    <- "black"  
axis_line_width    <- 0.8      
axis_tick_color    <- "black"  
axis_tick_width    <- 0.8      

label_host <- "Host Background"
label_mgs  <- "Host MGS"       
label_phag <- "Active Phage"

color_host_bg      <- "grey60"
color_mgs_bg       <- "#377EB8"
color_phage_act    <- "#E41A1C"

custom_colors      <- setNames(c(color_host_bg, color_mgs_bg, color_phage_act), 
                               c(label_host, label_mgs, label_phag))

m1_bg_color        <- "#E0ECF4"  
m1_line_width      <- 0.7        

m2_density_alpha   <- 0.5       
m2_density_line_col<- "black"    
m2_density_line_w  <- 0.5        

if (!dir.exists(out_dir_combined)) dir.create(out_dir_combined, recursive = TRUE)


get_interval <- function(v) {
  if (is.na(v)) return(NA_character_)
  if (v >= 0 & v < 0.47) return("0_0.47")
  if (v >= 0.47 & v < 0.70) return("0.47_0.70")
  if (v >= 0.70 & v <= 1.0) return("0.70_1.0")
  return(NA_character_)
}

load_and_merge_spa <- function(sample_prefix, source_type) {
  p_file <- file.path(csv_dir, paste0(sample_prefix, "_spa_precise.csv"))
  s_file <- file.path(csv_dir, paste0(sample_prefix, "_spa_sensitive.csv"))
  if (!file.exists(p_file) | !file.exists(s_file)) return(NULL)
  
  df_p <- suppressWarnings(fread(p_file)) %>% mutate(Delta_Mean = as.numeric(Delta_Mean))
  df_s <- suppressWarnings(fread(s_file)) %>% mutate(Effect_Size_Value = as.numeric(Effect_Size_Value))
    
  inner_join(df_p %>% select(Folder, Genome, Phage_File, Delta_Mean), 
             df_s %>% select(Folder, Genome, Phage_File, Effect_Size_Value), 
             by = c("Folder", "Genome", "Phage_File")) %>%
    mutate(
      Int_P = sapply(Delta_Mean, get_interval),
      Int_S = sapply(Effect_Size_Value, get_interval),
      Source_Type = source_type
    ) %>%
    filter(!is.na(Int_P) & !is.na(Int_S) & Int_P == Int_S) %>%
    rename(Final_Interval = Int_P)
}

target_list <- bind_rows(
  load_and_merge_spa("Con_AB", "AB"),
  load_and_merge_spa("M_AB", "AB"),
  load_and_merge_spa("GBIC", "GBIC")
)

shared_theme <- theme_classic(base_family = global_base_family) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 12, color = "black"),
    axis.title.x = element_text(size = size_axis_title, color = color_axis_title, face = face_axis_title),
    axis.title.y = element_text(size = size_axis_title, color = color_axis_title, face = face_axis_title),
    axis.text.x  = element_text(size = size_axis_text, color = color_axis_text, face = face_axis_text),
    axis.text.y  = element_text(size = size_axis_text, color = color_axis_text, face = face_axis_text),
    axis.line    = element_line(color = axis_line_color, linewidth = axis_line_width),
    axis.ticks   = element_line(color = axis_tick_color, linewidth = axis_tick_width),
    legend.position = legend_position,
    legend.title    = element_blank(),
    legend.text     = element_text(size = legend_text_size, color = legend_text_color, face = legend_text_face),
    legend.background = element_blank()
  )

for (intv in c("0_0.47", "0.47_0.70", "0.70_1.0")) {
  dir.create(file.path(out_dir_combined, intv), recursive = TRUE, showWarnings = FALSE)
}


process_phage <- function(i) {
  row_data <- target_list[i, ]
  folder   <- row_data$Folder
  genome   <- row_data$Genome
  phage_f  <- row_data$Phage_File
  src_type <- row_data$Source_Type
  interval <- row_data$Final_Interval
  
  clean_phage_name <- gsub("phage_depth_|\\.txt", "", phage_f)
  
  p1 <- ggplot() + theme_void() + ggtitle("Profile Data Missing") + theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggplot() + theme_void() + ggtitle("Density Data Missing") + theme(plot.title = element_text(hjust = 0.5))

  base_dir <- if (src_type == "AB") ab_depth_base else gbic_depth_base
  data_dir <- file.path(base_dir, folder, "depth_to_stat", genome)
  
  host_file  <- file.path(data_dir, paste0(genome, "_host_nonphage_depth.txt"))
  phage_file <- file.path(data_dir, phage_f)
  mgs_file   <- file.path(data_dir, paste0(genome, "_mgs_depth.txt"))
  
  if (file.exists(host_file) && file.exists(phage_file)) {
    
    dt_host_full <- suppressWarnings(fread(host_file, header = FALSE, col.names = c("contig", "position", "depth")))
    dt_phage     <- suppressWarnings(fread(phage_file, header = FALSE, col.names = c("contig", "position", "depth")))
    

    if (file.exists(mgs_file)) {
      dt_mgs <- suppressWarnings(fread(mgs_file, header = FALSE, col.names = c("contig", "position", "depth")))
    } else {
      dt_mgs <- data.table(depth = numeric(), Group = character())
    }
    
    if (nrow(dt_host_full) > 0 && nrow(dt_phage) > 0) {
   
      dt_host_full[, Group := label_host]
      dt_phage[, Group := label_phag]
      
      dt_all <- rbind(dt_host_full, dt_phage)
      
      dt_all[, contig_number := as.numeric(str_extract(contig, "\\d+$"))]
      if (any(is.na(dt_all$contig_number))) {
        dt_all[, contig_number := as.numeric(as.factor(contig))]
      }
      setorder(dt_all, contig_number, position)
      
      dt_all[, global_row := .I]
      
      target_contig <- unique(dt_phage$contig)[1]
      phage_start_pos <- min(dt_phage$position)
      phage_end_pos   <- max(dt_phage$position)
      
      start_row <- dt_all[contig == target_contig & position == phage_start_pos, global_row][1]
      end_row   <- dt_all[contig == target_contig & position == phage_end_pos, global_row][1]
      
      if (!is.na(start_row) && !is.na(end_row)) {
        
        buffer_size <- 30000
        plot_start_row <- max(1, start_row - buffer_size)
        plot_end_row   <- min(nrow(dt_all), end_row + buffer_size)
        
        dt_profile <- dt_all[global_row >= plot_start_row & global_row <= plot_end_row]
        
        dt_profile[, line_group := case_when(
          Group == label_phag ~ "phage",
          global_row < start_row ~ "host_left",
          global_row > end_row ~ "host_right",
          TRUE ~ "host_mid"
        )]
        
        p1 <- ggplot(dt_profile, aes(x = global_row, y = depth)) +
          annotate("rect", xmin = start_row, xmax = end_row, ymin = 0, ymax = Inf, fill = m1_bg_color, alpha = 1) +
          geom_line(aes(color = Group, group = line_group), linewidth = m1_line_width, alpha = 0.8) +
          scale_color_manual(values = custom_colors) +
          labs(title = "A. Genomic Depth Profile (Global Index)", 
               x = "Global Position Index (Concatenated Contigs)", 
               y = "Depth") + 
          shared_theme
      }
      

      if (is.na(var(dt_phage$depth)) || var(dt_phage$depth) == 0) {
        dt_phage$depth[1] <- dt_phage$depth[1] + 1e-6; dt_phage$depth[2] <- dt_phage$depth[2] - 1e-6
      }
      
      plot_data2 <- dt_all[, .(depth, Group)] 
      
 
      if (nrow(dt_mgs) > 0) {
        dt_mgs[, Group := label_mgs]
        plot_data2 <- rbind(plot_data2, dt_mgs[, .(depth, Group)])
      }
      
    
      plot_data2$Group <- factor(plot_data2$Group, levels = c(label_host, label_mgs, label_phag))
      
      plot_data2$depth <- plot_data2$depth + 1
      
      p2 <- ggplot(plot_data2, aes(x = depth, fill = Group)) +
        geom_density(color = m2_density_line_col, alpha = m2_density_alpha, linewidth = m2_density_line_w) +
        scale_x_log10(labels = scales::label_scientific()) +
        scale_fill_manual(values = custom_colors) +
        labs(title = "B. Depth Density Distribution (Log10 scale)", x = "Depth (count+1)", y = "Density") + 
        shared_theme
    }
  }
  

  combined_plot <- (p1 | p2) +
    plot_annotation(
      title = paste0("Phage Depth Overview: ", folder, " | ", clean_phage_name),
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  out_name <- paste0(folder, "_", clean_phage_name, "_dashboard.pdf")

  ggsave(file.path(out_dir_combined, interval, out_name), combined_plot, width = plot_combined_w, height = plot_combined_h, dpi = global_dpi, device = cairo_pdf)
  
  return(out_name)
}

results <- mclapply(1:nrow(target_list), process_phage, mc.cores = num_threads)

