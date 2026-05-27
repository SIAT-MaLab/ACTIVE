suppressPackageStartupMessages({
  library(tidyverse)
  library(tidygraph)
  library(ggraph)
  library(scales)
})

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

NODE_STROKE_WIDTH_SIG <- 0.8
NODE_STROKE_WIDTH_REG <- 0.2

PIE_BORDER_COLOR <- "black"

NET_HOST_SIZE_RANGE <- c(3, 15)
VIRUS_NODE_SIZE     <- 2.5

NET_EDGE_COLOR <- "grey40"
NET_EDGE_WIDTH <- 0.4
NET_EDGE_ALPHA <- 0.4

LABEL_TEXT_SIZE  <- 3.5
LABEL_TEXT_COLOR <- "black"
LABEL_FONT_FACE  <- "bold"

PIE_WIDTH  <- 8;  PIE_HEIGHT <- 8
NET_WIDTH  <- 18; NET_HEIGHT <- 14

message(">>> [1/5] Loading full metadata and Exact test statistical results...")

full_meta <- read_tsv("merged_phage_stats_taxonomy.tsv", show_col_types = FALSE)
stats_exact <- read_tsv("stats_baseline_vOTU_EXACT.tsv", show_col_types = FALSE)

full_meta <- full_meta %>%
  mutate(Activity_Score = as.numeric(as.character(Activity_Score)))

active_votu_list <- full_meta %>%
  filter(Activity_Score >= 0.7) %>%
  filter(!is.na(vOTU) & vOTU != "") %>%
  pull(vOTU) %>%
  unique()

message(sprintf(">>> Identified %d Active vOTUs.", length(active_votu_list)))

message(">>> [2/5] Recalculating true cross-host categories for all active vOTUs...")

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

votu_sig_map <- stats_exact %>%
  select(vOTU, Significance_Type, p.adj) %>%
  distinct() %>%
  mutate(Is_Significant = (Significance_Type == "Significant"))

active_votu_info <- tibble(vOTU = active_votu_list) %>%
  left_join(votu_cross_full, by = "vOTU") %>%
  left_join(votu_sig_map, by = "vOTU") %>%
  mutate(
    Cross_Category = replace_na(Cross_Category, "Non-crossing"),
    Is_Significant = replace_na(Is_Significant, FALSE)
  )

message(">>> [3/5] Generating Cross-Host Potential Pie Chart...")

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

message(">>> [4/5] Extracting edges and constructing Host-Phage Network (Host = Family)...")

edges <- full_meta %>%
  filter(vOTU %in% active_votu_list) %>%
  filter(!is.na(Family) & Family != "None" & Family != "") %>%
  select(vOTU, Host = Family) %>%
  distinct()

host_colors_map <- full_meta %>%
  filter(!is.na(Family) & Family != "None" & Family != "") %>%
  select(Family, Phylum) %>%
  distinct(Family, .keep_all = TRUE)

graph <- as_tbl_graph(edges, directed = FALSE) %>%
  activate(nodes) %>%
  mutate(
    node_type = ifelse(name %in% active_votu_list, "Virus", "Host"),
    degree = centrality_degree()
  ) %>%
  left_join(host_colors_map, by = c("name" = "Family")) %>%
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
    host_rank = ifelse(node_type == "Host", rank(-degree, ties.method = "first", na.last = "keep"), NA),
    virus_rank = ifelse(node_type == "Virus" & Is_Significant == TRUE, rank(p.adj, ties.method = "first", na.last = "keep"), NA),
    show_label = case_when(
      node_type == "Host" & host_rank <= 10 ~ TRUE,
      node_type == "Virus" & virus_rank <= 10 ~ TRUE,
      TRUE ~ FALSE
    ),
    node_label = ifelse(show_label, name, "")
  )

fill_palette <- c(PHYLUM_COLORS, CROSS_HOST_COLORS)
missing_groups <- setdiff(unique(graph %>% as_tibble() %>% pull(color_group)), names(fill_palette))
if(length(missing_groups) > 0) {
  fill_palette <- c(fill_palette, setNames(rep("#666666", length(missing_groups)), missing_groups))
}

border_palette <- c(
  "Host"         = BORDER_HOST,
  "Virus_Sig"    = BORDER_VIRUS_SIG,
  "Virus_NonSig" = BORDER_VIRUS_NONSIG
)

message(">>> [5/5] Rendering Network Plot (Applying Stress layout algorithm, please wait)...")

p_net <- ggraph(graph, layout = 'stress') +
  geom_edge_link(color = NET_EDGE_COLOR, width = NET_EDGE_WIDTH, alpha = NET_EDGE_ALPHA) +
  geom_node_point(
    data = function(x) filter(x, node_type == "Host"),
    aes(fill = color_group, size = degree, shape = node_type,
        color = border_group, stroke = stroke_width)
  ) +
  geom_node_point(
    data = function(x) filter(x, node_type == "Virus"),
    aes(fill = color_group, shape = node_type,
        color = border_group, stroke = stroke_width),
    size = VIRUS_NODE_SIZE
  ) +
  geom_node_text(aes(label = node_label),
                 repel = TRUE,
                 size = LABEL_TEXT_SIZE,
                 color = LABEL_TEXT_COLOR,
                 fontface = LABEL_FONT_FACE,
                 min.segment.length = 0,
                 segment.color = "grey40",
                 segment.size = 0.4,
                 max.overlaps = Inf) +
  scale_shape_manual(name = "Node Type", values = c("Host" = 22, "Virus" = 21)) +
  scale_size_continuous(range = NET_HOST_SIZE_RANGE, name = "Host Degree") +
  scale_fill_manual(values = fill_palette, name = "Phylum / True Cross-Category") +
  scale_color_manual(values = border_palette, name = "Significance (Border)") +
  scale_discrete_identity(aesthetics = "stroke") +
  theme_graph() +
  labs(
    title = "True Potential Network of Active vOTUs (Host = Family)",
    subtitle = "Edges: Full host range. Viruses: Red border = Exact test SAvOTU. Labels: Top 10 Host Families & Top 10 SAvOTUs."
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
    color = guide_legend(override.aes = list(size = 6, shape = 21, fill = "white", stroke = 1.5))
  )

out_file_pie <- "Active_vOTU_Pie_CrossHost_Categories.pdf"
ggsave(out_file_pie, p_pie, width = PIE_WIDTH, height = PIE_HEIGHT, bg = "white", device = cairo_pdf)
message(">>> Pie chart saved successfully to: ", out_file_pie)

out_file_net <- "Active_vOTU_Network_Exact_Family_Level.pdf"
ggsave(out_file_net, p_net, width = NET_WIDTH, height = NET_HEIGHT, bg = "white", device = cairo_pdf)
message(">>> Network graph saved successfully to: ", out_file_net)

message(">>> All tasks completed successfully!")
