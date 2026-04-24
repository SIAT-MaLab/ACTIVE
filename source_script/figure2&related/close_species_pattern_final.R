
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggsci)
  library(ggpubr)
})

pca_file   <- "../species_pca_coordinates.tsv"
skani_file <- "../species_skani_distances.tsv"
out_img    <- "plot_output_merged/Result_Final_Scatter_Plot_Professional.pdf"

custom_phylum_colors <- c(
  "Actinomycetota"    = "#f5cfa6",
  "Pseudomonadota"    = "#e68f9f",  
  "Bacteroidota"      = "#efe3ef",
  "Bacillota"         = "#7abbce",
  "Bacillota_A"       = "#dce8ef",
  "Verrucomicrobiota" = "#d4d68a",
  "Bacillota_I"       = "#34666b",
  "Bacillota_C"       = "#add2e5",
  "Bacillota_B"       = "#769499",
  "Fusobacteriota"    = "#231815",
  "Desulfobacterota"  = "#d5d5d6"
)

point_size      <- 2.5       
point_alpha     <- 1      

reg_line_color  <- "black"  
reg_line_width  <- 1.5      
reg_ci_fill     <- "gray70"  
reg_ci_alpha    <- 0.6       

global_dpi         <- 300   
plot_width         <- 12     
plot_height        <- 8     
global_font_family <- 

panel_border_color <- "black"
panel_border_width <- 1.5    

axis_title_size    <- 20
axis_title_color   <- "black"
axis_title_face    <- "bold"  

axis_text_size     <- 18
axis_text_color    <- "black"
axis_text_face     <- "plain"

legend_title_size  <- 16
legend_title_face  <- "bold"
legend_text_size   <- 11
legend_text_face   <- "plain"
legend_position    <- "right" 

setDTthreads(parallel::detectCores(logical = FALSE))

if(!file.exists(pca_file)) stop(paste("missing", pca_file))
dt_pca <- fread(pca_file)
setkey(dt_pca, Species)

if(!file.exists(skani_file)) stop(paste("missing", skani_file))
dt_skani_raw <- fread(skani_file)

list_pca   <- unique(dt_pca$Species)
list_skani <- unique(c(dt_skani_raw$Species_A, dt_skani_raw$Species_B))

diff_skani_names <- setdiff(list_skani, list_pca)
dt_diff_skani <- data.table(Species = diff_skani_names)
dt_diff_skani[, Reason := "Excluded by PCA Quality Control (Genome Count < 3)"]
dt_diff_skani[, Source := "Present in Skani, Absent in PCA"]

diff_pca_names <- setdiff(list_pca, list_skani)
dt_diff_pca <- dt_pca[Species %in% diff_pca_names, .(Species, Genus, Total)]
dt_diff_pca[, Reason := "Singleton Species in Genus (No other species available for pairwise comparison)"]
dt_diff_pca[, Source := "Present in PCA, Absent in Skani"]

fwrite(dt_diff_skani, "Report_Diff_Skani_Only.tsv", sep = "\t")
fwrite(dt_diff_pca,   "Report_Diff_PCA_Only.tsv", sep = "\t")

dt_agg <- dt_skani_raw[, .(
  Genomic_Dist_Mean = mean(Genomic_Dist),
  ANI_Mean = mean(ANI),
  N_Pairs = .N
), by = .(Genus, Species_A, Species_B)]

dt_agg[dt_pca, on = .(Species_A = Species), `:=`(PC1_A = i.PC1, PC2_A = i.PC2, Phylum = i.Phylum)]
dt_agg[dt_pca, on = .(Species_B = Species), `:=`(PC1_B = i.PC1, PC2_B = i.PC2)]

dt_final <- dt_agg[!is.na(PC1_A) & !is.na(PC1_B) & !is.na(Phylum)]
dt_final[, Pheno_Dist := sqrt((PC1_A - PC1_B)^2 + (PC2_A - PC2_B)^2)]

fwrite(dt_final, "Result_Final_Merged_Distances.tsv", sep = "\t")

plot_df <- as_tibble(dt_final)

stats_data <- plot_df %>%
  group_by(Phylum) %>%
  summarise(
    n_obs = n(),
    sd_x = sd(Genomic_Dist_Mean),
    sd_y = sd(Pheno_Dist),
    r_val = if(n() >= 3 && sd(Genomic_Dist_Mean) > 0 && sd(Pheno_Dist) > 0) 
              cor(Genomic_Dist_Mean, Pheno_Dist, method = "pearson") else NA_real_,
    p_val = if(n() >= 3 && sd(Genomic_Dist_Mean) > 0 && sd(Pheno_Dist) > 0) 
              cor.test(Genomic_Dist_Mean, Pheno_Dist)$p.value else NA_real_,
    .groups = 'drop'
  ) %>%
  filter(!is.na(r_val) & !is.na(p_val)) %>%
  mutate(
    stat_label = paste0(
      Phylum, ": R=", sprintf("%.2f", r_val), 
      " P=", ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val))
    )
  )

print(stats_data$stat_label)

unique_phyla <- unique(plot_df$Phylum)
my_palette <- c()
missing_phyla <- c()

for (p in unique_phyla) {
  if (p %in% names(custom_phylum_colors)) {
    my_palette[p] <- custom_phylum_colors[p]
  } else {
    missing_phyla <- c(missing_phyla, p)
  }
}

if (length(missing_phyla) > 0) {
  base_cols <- pal_npg("nrc")(min(length(missing_phyla), 10))
  fallback_cols <- if(length(missing_phyla) > length(base_cols)) colorRampPalette(base_cols)(length(missing_phyla)) else base_cols[1:length(missing_phyla)]
  for (i in seq_along(missing_phyla)) {
    my_palette[missing_phyla[i]] <- fallback_cols[i]
  }
}

p_final <- ggplot(plot_df, aes(x = Genomic_Dist_Mean, y = Pheno_Dist)) +

  geom_point(aes(color = Phylum), alpha = point_alpha, size = point_size) +
  scale_color_manual(values = my_palette) + 

  geom_smooth(aes(group = 1), method = "lm", 
              color = reg_line_color, linewidth = reg_line_width, 
              fill = reg_ci_fill, alpha = reg_ci_alpha) +

  geom_point(data = stats_data, 
             aes(x = min(plot_df$Genomic_Dist_Mean, na.rm = TRUE), 
                 y = min(plot_df$Pheno_Dist, na.rm = TRUE), 
                 shape = stat_label), 
             alpha = 0) +

  guides(
    color = guide_legend(title = "Phylum Group", order = 1, override.aes = list(size = 4, alpha = 1)),
    shape = guide_legend(title = "Statistics (Valid Groups)", order = 2, override.aes = list(alpha = 1, size = 0, stroke = 0)) 
  ) +

  labs(
    title = "Evolutionary Divergence Landscape",
    subtitle = paste0("Total Species Pairs: ", nrow(plot_df), " (Overall trend in black)"),
    x = "Genomic Distance (1 - Average ANI)",
    y = "Lysogeny Phenotypic Distance (PCA Euclidean)"
  ) +

  theme_bw(base_family = global_font_family) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = panel_border_color, linewidth = panel_border_width, fill = NA),
    axis.title.x = element_text(size = axis_title_size, color = axis_title_color, face = axis_title_face, margin = margin(t = 10)),
    axis.title.y = element_text(size = axis_title_size, color = axis_title_color, face = axis_title_face, margin = margin(r = 10)),
    axis.text.x = element_text(size = axis_text_size, color = axis_text_color, face = axis_text_face),
    axis.text.y = element_text(size = axis_text_size, color = axis_text_color, face = axis_text_face),
    legend.position = legend_position,
    legend.box = "vertical",
    legend.background = element_blank(),
    legend.title = element_text(size = legend_title_size, face = legend_title_face, color = "black"),
    legend.text = element_text(size = legend_text_size, face = legend_text_face, color = "black"),
    plot.title = element_text(size = axis_title_size + 2, face = "bold", color = "black"),
    plot.subtitle = element_text(size = axis_title_size - 1, color = "gray30")
  )

ggsave(out_img, p_final, width = plot_width, height = plot_height, dpi = global_dpi, device = cairo_pdf)

