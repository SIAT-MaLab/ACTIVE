
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggsci)
  library(patchwork)
})

input_file       <- "../prophage_stats_final.tsv"
out_dir          <- "./plot_output_merged"

min_host_count   <- 5      
threshold_count  <- 100    
threshold_active <- 0.20   

target_levels    <- c("Phylum", "Species") 

custom_phylum_colors <- c(
  "Actinomycetota"    = "#61CBA6",  
  "Pseudomonadota"    = "#F8767D",  
  "Bacteroidota"      = "#6166F8",  
  "Bacillota"         = "#D49A00",
  "Bacillota_A"       = "#C66500",
  "Verrucomicrobiota" = "#7425C2",
  "Bacillota_I"       = "#F4DE61",
  "Bacillota_C"       = "#635100",
  "Bacillota_B"       = "#FFFF42",
  "Fusobacteriota"    = "#FF00D4",
  "Desulfobacterota"  = "#4CFFFF"
)

phylum_point_alpha <- 0.85
cols_bar <- c("Non" = "grey60", "Single" = "#E41A1C", "Poly" = "#74004A")
shape_levels <- c("Actinomycetota", "Bacillota*", "Bacteroidota", "Pseudomonadota", "Others")
custom_shapes <- c(
  "Actinomycetota" = 21, 
  "Bacillota*"     = 22, 
  "Bacteroidota"   = 23, 
  "Pseudomonadota" = 24, 
  "Others"         = 25 
)

color_other_bg <- "grey85" 
global_dpi        <- 600      
global_base_size  <- 20       
plot_phylum_w     <- 22
plot_phylum_h_min <- 10       

plot_species_w    <- 25
plot_species_h    <- 12

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

raw_data <- read_tsv(input_file, show_col_types = FALSE)
id_col_name <- colnames(raw_data)[1]

clean_data <- raw_data %>%
  mutate(Real_Genome_ID = str_remove(.data[[id_col_name]], "_[0-9]+$")) %>%
  group_by(Real_Genome_ID, Domain, Phylum, Class, Order, Family, Genus, Species, source) %>%
  summarise(
    Prophage_Count = sum(Prophage_Count, na.rm = TRUE),
    Active_Prophage_Count = sum(Active_Prophage_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Lysogeny_Status = case_when(
    Prophage_Count == 0 ~ "Non-lysogen",
    Prophage_Count == 1 ~ "Single-lysogen",
    Prophage_Count > 1  ~ "Poly-lysogen"
  ))

for (plot_level in target_levels) {
 
  df_subset <- clean_data %>% drop_na(all_of(plot_level))

  parent_col <- switch(plot_level, "Phylum" = "Domain", "Species" = "Genus")
  
  stats_df <- df_subset %>%
    group_by(Taxon = .data[[plot_level]], Parent = .data[[parent_col]], Phylum) %>%
    summarise(
      Total_Strains = n(),
      Count_Non = sum(Lysogeny_Status == "Non-lysogen"),
      Count_Single = sum(Lysogeny_Status == "Single-lysogen"),
      Count_Poly = sum(Lysogeny_Status == "Poly-lysogen"),
      Count_Lysogens = sum(Prophage_Count > 0),
      Count_Active_Hosts = sum(Active_Prophage_Count > 0),
      .groups = "drop"
    ) %>%
    filter(Total_Strains >= min_host_count) %>%
    mutate(
      Rate_Lysogeny_X = Count_Lysogens / Total_Strains,
      Rate_Active_Host_Y = ifelse(Count_Lysogens > 0, Count_Active_Hosts / Count_Lysogens, 0),
      Ratio_Non = Count_Non / Total_Strains, 
      Ratio_Single = Count_Single / Total_Strains, 
      Ratio_Poly = Count_Poly / Total_Strains
    )
  
  if (nrow(stats_df) == 0) next

  global_mean_lysogeny <- sum(stats_df$Count_Lysogens) / sum(stats_df$Total_Strains)
  total_lysogens_pop <- sum(stats_df$Count_Lysogens)
  global_mean_active <- if(total_lysogens_pop > 0) sum(stats_df$Count_Active_Hosts) / total_lysogens_pop else 0

  stats_df <- stats_df %>%
    mutate(is_highlight = Total_Strains >= threshold_count | Rate_Active_Host_Y >= threshold_active)

  if (plot_level == "Phylum") {
    stats_df <- stats_df %>% mutate(Color_Group = ifelse(is_highlight, Taxon, "Other"))
    legend_title_text <- "Phylum"
  } else {
    highlight_parents <- unique(stats_df$Parent[stats_df$is_highlight])
    stats_df <- stats_df %>% mutate(Color_Group = ifelse(Parent %in% highlight_parents, Parent, "Other"))
    legend_title_text <- paste0("Taxonomy (", parent_col, ")")
  }
  
  sorted_groups <- sort(setdiff(unique(stats_df$Color_Group), "Other"))
  stats_df$Color_Group <- factor(stats_df$Color_Group, levels = c("Other", sorted_groups))
  stats_df <- stats_df %>% arrange(Color_Group)

  my_palette <- c("Other" = color_other_bg)
  n_needed <- length(sorted_groups)
  
  if (n_needed > 0) {
    if (plot_level == "Phylum") {
      predefined <- sorted_groups[sorted_groups %in% names(custom_phylum_colors)]
      missing <- sorted_groups[!sorted_groups %in% names(custom_phylum_colors)]
      for (p in predefined) my_palette[p] <- custom_phylum_colors[p]
      if (length(missing) > 0) {
        base_cols <- pal_npg("nrc")(min(length(missing), 10))
        final_cols <- if(length(missing) > length(base_cols)) colorRampPalette(base_cols)(length(missing)) else base_cols[1:length(missing)]
        for (i in seq_along(missing)) my_palette[missing[i]] <- final_cols[i]
      }
    } else {

      base_cols <- pal_npg("nrc")(min(n_needed, 10))
      final_cols <- if(n_needed > length(base_cols)) colorRampPalette(base_cols)(n_needed) else base_cols[1:n_needed]
      my_palette <- c("Other" = color_other_bg, setNames(final_cols, sorted_groups))
    }
  }
  
  export_df <- stats_df %>%
    mutate(
      Global_Mean_X = global_mean_lysogeny, 
      Global_Mean_Y = global_mean_active,
      Quadrant = case_when(
        Rate_Lysogeny_X >= global_mean_lysogeny & Rate_Active_Host_Y >= global_mean_active ~ "Q1",
        Rate_Lysogeny_X <  global_mean_lysogeny & Rate_Active_Host_Y >= global_mean_active ~ "Q2",
        Rate_Lysogeny_X <  global_mean_lysogeny & Rate_Active_Host_Y <  global_mean_active ~ "Q3",
        Rate_Lysogeny_X >= global_mean_lysogeny & Rate_Active_Host_Y <  global_mean_active ~ "Q4"
      )
    )
  
  csv_out <- file.path(out_dir, paste0("Statistics_", plot_level, ".csv"))
  write_csv(export_df, csv_out)

  if (plot_level == "Phylum") {
    
    stats_df <- stats_df %>% mutate(Label_Text = ifelse(is_highlight, Taxon, NA))
    p1 <- ggplot(stats_df, aes(x = Rate_Lysogeny_X, y = Rate_Active_Host_Y)) +
      geom_vline(xintercept = global_mean_lysogeny, linetype = "dashed", color = "black", alpha=0.6) +
      geom_hline(yintercept = global_mean_active, linetype = "dashed", color = "black", alpha=0.6) +
      geom_point(aes(size = Total_Strains, fill = Color_Group), 
                 shape = 21, color = "black", stroke = 0.3, alpha = phylum_point_alpha) +
      geom_text_repel(aes(label = Label_Text, color = Color_Group), 
                      size = 5, fontface = "bold", bg.color = "white", bg.r = 0.15,
                      box.padding = 0.6, max.overlaps = Inf, show.legend = FALSE) +
      scale_fill_manual(values = my_palette, name = legend_title_text) +
      scale_color_manual(values = my_palette, guide = "none") +
      scale_size_continuous(trans = "log10", range = c(3, 12), breaks = c(10, 100, 1000, 10000), 
                            labels = c("10", "100", "1k", "10k"), name = "Genome Count") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)), limits = c(0, max(stats_df$Rate_Active_Host_Y) * 1.05)) +
      scale_x_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0.01, 0.05))) +
      labs(title = paste0("A. Overview (", plot_level, ")"), 
           subtitle = "Quadrants defined by population weighted means",
           x = "Lysogen Rate", y = "Active Rate") +
      theme_bw(base_size = global_base_size) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        plot.title = element_text(face = "bold", size = 16),
        legend.title = element_text(hjust = 0.5) 
      ) +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, order = 1),
             size = guide_legend(title.position = "top", title.hjust = 0.5, order = 2))
    df_hi <- stats_df %>% filter(is_highlight)
    df_oth <- stats_df %>% filter(!is_highlight)
    if(nrow(df_oth) > 0) df_hi <- bind_rows(df_hi, df_oth %>% summarise(Taxon="Other", Total_Strains=sum(Total_Strains), Ratio_Non=0, Ratio_Single=0, Ratio_Poly=0))
    
    levs <- df_hi %>% filter(Taxon!="Other") %>% arrange(Ratio_Non) %>% pull(Taxon)
    if("Other" %in% df_hi$Taxon) levs <- c("Other", levs)
    df_hi$Taxon <- factor(df_hi$Taxon, levels = levs)
    
    bar_long <- df_hi %>% select(Taxon, Ratio_Non, Ratio_Single, Ratio_Poly) %>% 
      pivot_longer(-Taxon, names_to="Type", values_to="Val") %>%
      mutate(Type = factor(str_remove(Type,"Ratio_"), levels=c("Non","Single","Poly")))
    p2 <- ggplot(bar_long, aes(Taxon, Val, fill=Type)) + 
      geom_col(color="black", width=0.8, linewidth=0.2) + coord_flip() +
      scale_fill_manual(values=cols_bar) + 
      scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
      theme_bw(base_size=global_base_size) + 
      labs(title="B. Lysogeny Patterns", x=NULL, y="Proportion") +
      theme(legend.position="bottom", axis.text.y = element_text(face="italic"), legend.title = element_text(hjust = 0.5)) +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
    p3 <- ggplot(df_hi, aes(Taxon, Total_Strains)) + 
      geom_col(fill="grey80", color="black", width=0.8, linewidth=0.2) + 
      geom_text(aes(label=Total_Strains), hjust=-0.1, size=4.5) + coord_flip() + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + 
      theme_bw(base_size=global_base_size) +
      labs(title="C. Count", x=NULL, y="Genomes") + 
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid=element_blank())
    combined_plot <- p1 + p2 + p3 + plot_layout(widths = c(1.5, 1, 0.4))
    plot_h <- max(plot_phylum_h_min, 0.4 * nrow(df_hi) + 4) 
    png_out <- file.path(out_dir, paste0("Merged_Analysis_Quadrants_", plot_level, ".png"))
    ggsave(png_out, combined_plot, width = plot_phylum_w, height = plot_h, dpi = global_dpi, limitsize = FALSE)
    message(sprintf("  -> 已保存图片: %s", png_out))
    
  } else if (plot_level == "Species") {
    stats_df <- stats_df %>%
      mutate(Shape_Group = case_when(
        Phylum == "Actinomycetota" ~ "Actinomycetota",
        str_detect(Phylum, "^Bacillota") ~ "Bacillota*",
        Phylum == "Bacteroidota" ~ "Bacteroidota",
        Phylum == "Pseudomonadota" ~ "Pseudomonadota",
        TRUE ~ "Others"
      ))
    stats_df$Shape_Group <- factor(stats_df$Shape_Group, levels = shape_levels)
    
    p1 <- ggplot(stats_df, aes(x = Rate_Lysogeny_X, y = Rate_Active_Host_Y)) +
      geom_vline(xintercept = global_mean_lysogeny, linetype = "dashed", color = "black", alpha=0.6) +
      geom_hline(yintercept = global_mean_active, linetype = "dashed", color = "black", alpha=0.6) +
      geom_point(aes(size = Total_Strains, fill = Color_Group, shape = Shape_Group), 
                 color = "black", stroke = 0.3, alpha = phylum_point_alpha) +
      scale_shape_manual(values = custom_shapes, name = "Phylum Group") +
      scale_fill_manual(values = my_palette, name = legend_title_text) +
      scale_size_continuous(trans = "log10", range = c(3, 12), breaks = c(10, 100, 1000, 10000), 
                            labels = c("10", "100", "1k", "10k"), name = "Genome Count") +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)), limits = c(0, max(stats_df$Rate_Active_Host_Y) * 1.05)) +
      scale_x_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0.01, 0.05))) +
      labs(title = paste0("Overview (", plot_level, ")"), 
           subtitle = paste0("Mean: X=", round(global_mean_lysogeny,2), ", Y=", round(global_mean_active,2)),
           x = "Lysogen Rate", y = "Active Rate") +
      theme_bw(base_size = 16)
    
    legend_style <- theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.justification = "top",
      legend.key.size = unit(1.2, "cm"), 
      legend.text = element_text(size = 14), 
      legend.title = element_text(size = 15, face = "bold"),
      legend.spacing.y = unit(0.5, "cm")
    )
    
    n_taxa <- length(unique(stats_df$Color_Group))
    col_ncol <- if (n_taxa > 30) 3 else 1 
    
    guides_style <- guides(
      fill = guide_legend(title = legend_title_text, ncol = col_ncol, override.aes = list(size = 8, shape = 21), order = 1),
      shape = guide_legend(title = "Phylum Group", ncol = 1, override.aes = list(size = 8, fill = "grey60"), order = 2),
      size = guide_legend(title = "Count", ncol = 1, order = 3)
    )
    
    p1_final <- p1 + legend_style + guides_style
    
    png_out <- file.path(out_dir, paste0("Merged_Analysis_", plot_level, ".png"))
    ggsave(png_out, p1_final, width = plot_species_w, height = plot_species_h, dpi = global_dpi, limitsize = FALSE)

  }
}

