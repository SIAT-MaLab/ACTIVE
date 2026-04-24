
library(tidyverse)
library(ggrepel)
library(patchwork)

input_file <- "../prophage_stats_final.tsv"
min_genomes_per_species <- 3

target_phyla <- c("Actinomycetota", "Pseudomonadota", "Bacteroidota")

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
  "Fusobacteriota"    = "#231815"
)

raw_data <- read_tsv(input_file, show_col_types = FALSE)
id_col_name <- colnames(raw_data)[1]

clean_data_base <- raw_data %>%
  mutate(Real_Genome_ID = str_remove(.data[[id_col_name]], "_[0-9]+$"))

aggr_data <- clean_data_base %>%
  filter(!is.na(Species), Species != "NA", Species != "s__") %>%
  group_by(Real_Genome_ID, Family, Genus, Species, Phylum) %>%
  summarise(P_Count=sum(Prophage_Count,na.rm=T), A_Count=sum(Active_Prophage_Count,na.rm=T), .groups="drop")

features <- aggr_data %>%
  group_by(Species, Genus, Phylum, Family) %>%
  summarise(
    Total=n(),
    Rate_Lys=sum(P_Count>0)/n(),
    Rate_Act=ifelse(sum(P_Count>0)>0, sum(A_Count>0)/sum(P_Count>0), 0),
    Rate_Poly=sum(P_Count>1)/n(),
    Mean_Count=mean(P_Count),
    .groups="drop"
  ) %>%
  filter(Total >= min_genomes_per_species) %>%
  ungroup() %>%
  distinct(Rate_Lys, Rate_Act, Rate_Poly, Mean_Count, .keep_all = TRUE)

pca_mat <- features %>% select(Rate_Lys, Rate_Act, Rate_Poly, Mean_Count) %>% as.data.frame()
rownames(pca_mat) <- features$Species
pca_res <- prcomp(pca_mat, scale. = TRUE)
var_exp <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)

pca_coords <- as.data.frame(pca_res$x)
pca_coords$Species <- rownames(pca_coords)
plot_data <- left_join(pca_coords, features, by = "Species")

output_coord_file <- "species_pca_coordinates.tsv"
export_df <- plot_data %>%
  select(Species, Genus, Family, Phylum, PC1, PC2,
         Rate_Lys, Rate_Act, Rate_Poly, Mean_Count, Total)

write_tsv(export_df, output_coord_file)

loadings <- as.data.frame(pca_res$rotation)
loadings$Feature <- rownames(loadings)
arrow_scale <- 2.5
loadings$v1 <- loadings$PC1 * arrow_scale
loadings$v2 <- loadings$PC2 * arrow_scale

plot_data <- plot_data %>%
  mutate(Shape_Group = case_when(
    Phylum == "Actinomycetota" ~ "Actinomycetota",
    str_detect(Phylum, "^Bacillota") ~ "Bacillota*",
    Phylum == "Bacteroidota" ~ "Bacteroidota",
    Phylum == "Pseudomonadota" ~ "Pseudomonadota",
    TRUE ~ "Others"
  ))
shape_levels <- c("Actinomycetota", "Bacillota*", "Bacteroidota", "Pseudomonadota", "Others")
plot_data$Shape_Group <- factor(plot_data$Shape_Group, levels = shape_levels)
custom_shapes <- c(21, 22, 23, 24, 25) 

valid_phyla_all <- plot_data %>% count(Phylum) %>% filter(n >= 3) %>% pull(Phylum)
ellipse_data_all <- plot_data %>% filter(Phylum %in% valid_phyla_all)

valid_target <- plot_data %>% filter(Phylum %in% target_phyla) %>%
                count(Phylum) %>% filter(n >= 3) %>% pull(Phylum)
ellipse_data_selected <- plot_data %>% filter(Phylum %in% valid_target)

plot_pca_species <- function(ellipse_df, title_suffix) {
  p_scores <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey90") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey90") +
    stat_ellipse(data = ellipse_df, aes(color = Phylum),
                 type = "t", level = 0.90, linewidth = 2, show.legend = FALSE) +
    geom_point(aes(fill = Phylum, shape = Shape_Group, size = Total),
               color = "black", stroke = 0.3, alpha = 1) +
    scale_shape_manual(values = custom_shapes, name = "Phylum Group (Shape)") +
    scale_fill_manual(values = custom_phylum_colors, na.value = "grey50", name = "Phylum") +
    scale_color_manual(values = custom_phylum_colors, na.value = "grey50") +
    scale_size_continuous(range = c(1, 12), name = "Genome Count") +
    coord_fixed(ratio = 1) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 4), order = 1),
      shape = guide_legend(override.aes = list(size = 4, fill = "grey50"), order = 2),
      size = guide_legend(order = 3)
    ) +
    labs(title = paste0("A. Species Clustering (", title_suffix, ")"),
         x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw(base_size = 14)

  p_loadings <- ggplot() +
    annotate("path", x = cos(seq(0, 2*pi, len=100))*arrow_scale, y = sin(seq(0, 2*pi, len=100))*arrow_scale, color="grey80", linetype="dotted") +
    geom_vline(xintercept=0, linetype="dashed", color="grey90") + geom_hline(yintercept=0, linetype="dashed", color="grey90") +
    geom_segment(data = loadings, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.3,"cm")), color="#D32F2F", linewidth=1.2) +
    geom_text_repel(data = loadings, aes(x=v1, y=v2, label=Feature), color="black", fontface="bold", size=5, box.padding=0.5) +
    coord_fixed(ratio = 1) +
    labs(title = "B. Features", x="PC1", y="PC2") +
    theme_bw(base_size = 14) + theme(panel.grid = element_blank())

  return(p_scores + p_loadings + plot_layout(widths = c(1.8, 1)))
}

ggsave("Analysis_Species_PCA_AllPhyla.pdf", plot_pca_species(ellipse_data_all, "All Phyla"), width = 18, height = 8)
ggsave("Analysis_Species_PCA_Selected.pdf", plot_pca_species(ellipse_data_selected, "Selected Phyla"), width = 18, height = 8)

