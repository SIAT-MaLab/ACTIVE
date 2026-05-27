#!/usr/bin/env Rscript

library(tidyverse)
library(tidygraph)
library(ggraph)
library(scales)

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

CROSS_HOST_COLORS <- c(
  "Non-crossing"  = "grey80",
  "Cross-Species" = "#C0EBFF",
  "Cross-Genus"   = "#709DC5",
  "Cross-Family"  = "#E77C8E",
  "Cross-Order"   = "#D62828",
  "Cross-Class"   = "#9D0208",
  "Cross-Phylum"  = "#510014"
)

BORDER_VIRUS_SIG    <- "#861A0E" 
BORDER_VIRUS_NONSIG <- "grey80"  
BORDER_HOST         <- "black"   

NODE_STROKE_WIDTH_SIG <- 0.5
NODE_STROKE_WIDTH_REG <- 0.5

PIE_BORDER_COLOR <- "black"

NET_HOST_SIZE_RANGE <- c(3, 15)  
VIRUS_NODE_SIZE     <- 2.5       

NET_EDGE_COLOR <- "grey75"
NET_EDGE_WIDTH <- 0.4            
NET_EDGE_ALPHA <- 0.4            

LABEL_TEXT_SIZE  <- 3.5
LABEL_TEXT_COLOR <- "black"
LABEL_FONT_FACE  <- "bold"

PIE_WIDTH  <- 8;  PIE_HEIGHT <- 8
NET_WIDTH  <- 16; NET_HEIGHT <- 12 

full_meta <- read_tsv("../merged_phage_stats_taxonomy.tsv", show_col_types = FALSE)
stats_fisher <- read_tsv("../stats_fisher_vOTU.tsv", show_col_types = FALSE)


full_meta <- full_meta %>%
  mutate(Activity_Score = as.numeric(as.character(Activity_Score)))

active_votu_list <- full_meta %>%
  filter(Activity_Score >= 0.7) %>%
  filter(!is.na(vOTU) & vOTU != "") %>%
  pull(vOTU) %>%
  unique()

votu_cross_full <- full_meta %>%
  filter(vOTU %in% active_votu_list) %>%
  group_by(vOTU) %>%
  summarise(
    n_phylum  = n_distinct(Phylum, na.rm = TRUE),
    n_class   = n_distinct(Class, na.rm = TRUE),
    n_order   = n_distinct(Order, na.rm = TRUE),
    n_family  = n_distinct(Family, na.rm = TRUE),
    n_genus   = n_distinct(Genus, na.rm = TRUE),
    n_species = n_distinct(Species, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Cross_Category = case_when(
      n_phylum > 1  ~ "Cross-Phylum",
      n_class > 1   ~ "Cross-Class",
      n_order > 1   ~ "Cross-Order",
      n_family > 1  ~ "Cross-Family",
      n_genus > 1   ~ "Cross-Genus",
      n_species > 1 ~ "Cross-Species",
      TRUE          ~ "Non-crossing"
    )
  )

votu_sig_map <- stats_fisher %>%
  select(vOTU, p.adj) %>%
  distinct() %>%
  mutate(Is_Significant = (!is.na(p.adj) & p.adj <= 0.05))

active_votu_info <- tibble(vOTU = active_votu_list) %>%
  left_join(votu_cross_full, by = "vOTU") %>%
  left_join(votu_sig_map, by = "vOTU") %>%
  mutate(
    Cross_Category = replace_na(Cross_Category, "Non-crossing"),
    Is_Significant = replace_na(Is_Significant, FALSE)
  )

summary_pie <- active_votu_info %>%
  count(Cross_Category) %>%
  mutate(prop = n / sum(n))

p_pie <- ggplot(summary_pie, aes(x = "", y = n, fill = Cross_Category)) +
  geom_bar(stat = "identity", width = 1, color = PIE_BORDER_COLOR, linewidth = 0.5) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = CROSS_HOST_COLORS) +
  labs(
    title = "Proportion of True Cross-Host Potentials",
    subtitle = paste0("Calculated from Full Metadata for Active vOTUs (n = ", length(active_votu_list), ")")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

out_file_pie <- "Active_vOTU_Pie_TruePotential.pdf"
ggsave(out_file_pie, p_pie, width = PIE_WIDTH, height = PIE_HEIGHT, bg = "white", device = cairo_pdf)

edges <- full_meta %>%
  filter(vOTU %in% active_votu_list) %>%
  filter(!is.na(Genus) & Genus != "None" & Genus != "") %>%
  select(vOTU, Host = Genus) %>%
  distinct()

graph_global <- as_tbl_graph(edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(global_degree = centrality_degree()) %>%
  as_tibble()

host_global_degrees <- graph_global$global_degree[graph_global$name %in% edges$Host]
global_degree_limits <- c(min(host_global_degrees, na.rm = TRUE), max(host_global_degrees, na.rm = TRUE))

host_colors_map <- full_meta %>%
  filter(!is.na(Genus) & Genus != "None" & Genus != "") %>%
  select(Genus, Phylum) %>%
  distinct(Genus, .keep_all = TRUE)

stats_fisher_padj <- stats_fisher %>% select(vOTU, p.adj) %>% distinct()

host_colors_map <- host_colors_map %>%
  mutate(
    Plot_Group = case_when(
      Phylum %in% c("Pseudomonadota", "Bacteroidota") ~ "Pseudomonadota-Bacteroidota",
      grepl("^Bacillota", Phylum) ~ "Bacillota_Group", 
      Phylum == "Actinomycetota" ~ "Actinomycetota",
      TRUE ~ "Other"
    )
  )

unique_groups <- unique(host_colors_map$Plot_Group)

fill_palette <- c(PHYLUM_COLORS, CROSS_HOST_COLORS)
missing_groups <- setdiff(c(unique(host_colors_map$Phylum), unique(active_votu_info$Cross_Category)), names(fill_palette))
if(length(missing_groups) > 0) {
  fill_palette <- c(fill_palette, setNames(rep("#666666", length(missing_groups)), missing_groups))
}

border_palette <- c("Host" = BORDER_HOST, "Virus_Sig" = BORDER_VIRUS_SIG, "Virus_NonSig" = BORDER_VIRUS_NONSIG)


for(curr_group in unique_groups) {
  target_genera <- host_colors_map %>% filter(Plot_Group == curr_group) %>% pull(Genus)
  sub_edges <- edges %>% filter(Host %in% target_genera)
  if(nrow(sub_edges) == 0) next
  sub_graph <- as_tbl_graph(sub_edges, directed = FALSE) %>%
    activate(nodes) %>%
    mutate(node_type = ifelse(name %in% active_votu_list, "Virus", "Host")) %>%
    left_join(graph_global %>% select(name, global_degree), by = "name") %>%
    left_join(host_colors_map, by = c("name" = "Genus")) %>%
    left_join(active_votu_info, by = c("name" = "vOTU")) %>%
    mutate(
      color_group = case_when(
        node_type == "Host" ~ Phylum,
        node_type == "Virus" ~ as.character(Cross_Category)
      ),
      border_group = case_when(
        node_type == "Host" ~ "Host",
        node_type == "Virus" & Is_Significant == TRUE ~ "Virus_Sig",
        TRUE ~ "Virus_NonSig"
      ),
      stroke_width = ifelse(border_group == "Virus_Sig", NODE_STROKE_WIDTH_SIG, NODE_STROKE_WIDTH_REG)
    ) %>%
    mutate(
      host_rank = ifelse(node_type == "Host", rank(-global_degree, ties.method = "first", na.last = "keep"), NA),
      virus_rank = ifelse(node_type == "Virus" & Is_Significant == TRUE, rank(p.adj, ties.method = "first", na.last = "keep"), NA),
      
      show_label = case_when(
        node_type == "Host" & host_rank <= 10 ~ TRUE,
        node_type == "Virus" & virus_rank <= 10 ~ TRUE,
        TRUE ~ FALSE
      ),
      node_label = ifelse(show_label, name, "")
    )
  p_sub <- ggraph(sub_graph, layout = 'kk') +
    geom_edge_link(color = NET_EDGE_COLOR, width = NET_EDGE_WIDTH, alpha = NET_EDGE_ALPHA) +
    
    geom_node_point(data = function(x) filter(x, node_type == "Host"),
                    aes(fill = color_group, size = global_degree, shape = node_type, color = border_group, stroke = stroke_width)) +
    
    geom_node_point(data = function(x) filter(x, node_type == "Virus"),
                    aes(fill = color_group, shape = node_type, color = border_group, stroke = stroke_width),
                    size = VIRUS_NODE_SIZE) +
    
    geom_node_text(aes(label = node_label),
                   repel = TRUE, size = LABEL_TEXT_SIZE, color = LABEL_TEXT_COLOR, fontface = LABEL_FONT_FACE,
                   min.segment.length = 0, segment.color = "grey40", segment.size = 0.4, max.overlaps = Inf) +       
    
    scale_shape_manual(name = "Node Type", values = c("Host" = 22, "Virus" = 21)) +

    scale_size_continuous(range = NET_HOST_SIZE_RANGE, limits = global_degree_limits, name = "Global Host Degree") +
    scale_fill_manual(values = fill_palette, limits = names(fill_palette), drop = FALSE, name = "Phylum / True Cross-Category") +
    scale_color_manual(values = border_palette, limits = names(border_palette), drop = FALSE, name = "Significance (Border)") +
    scale_discrete_identity(aesthetics = "stroke") +
    
    theme_graph() +
    labs(
      title = paste("Network:", curr_group),
      subtitle = "Edges: Full host range. Viruses: Red border = Significant (p.adj <= 0.05). Labels: Local Top 10."
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 6, shape = 21, color = "black", stroke = 0.2)),
      shape = guide_legend(override.aes = list(size = 6, fill = "grey50", color = "black")),
      color = guide_legend(override.aes = list(size = 6, shape = 21, fill = "white", stroke = 1.5)),
      size = guide_legend()
    )
  

  out_file_name <- paste0("Network_Group_", curr_group, ".pdf")
  ggsave(out_file_name, p_sub, width = NET_WIDTH, height = NET_HEIGHT, bg = "white", device = cairo_pdf)
}

