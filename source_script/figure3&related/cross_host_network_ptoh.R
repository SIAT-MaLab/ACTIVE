library(tidyverse)
library(tidygraph)
library(ggraph)
library(ggsci)  
library(scales) 
library(forcats)

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
  "Fusobacteriota"    = "#b7a3c9"
)


data <- read_tsv("phage_full_stats_separated.tsv", show_col_types = FALSE)

votu_broad_status <- data %>%
  filter(!is.na(vOTU), !is.na(Species), Species != "None") %>%
  group_by(vOTU) %>%
  summarise(total_species_count = n_distinct(Species)) %>%
  mutate(is_broad = ifelse(total_species_count >= 10, "Broad", "Narrow"))

fam_to_phylum <- data %>%
  filter(!is.na(Family), !is.na(Phylum)) %>%
  select(Family, Phylum) %>%
  distinct(Family, .keep_all = TRUE)

genus_to_family <- data %>%
  filter(!is.na(Genus), !is.na(Family)) %>%
  select(Genus, Family) %>%
  distinct(Genus, .keep_all = TRUE)

host_taxonomy <- data %>%
  select(Species, Family, Phylum) %>%
  distinct(Species, .keep_all = TRUE)

get_phylum_palette <- function(groups) {
  host_groups <- setdiff(groups, c("Virus_Broad", "Virus_Shared", "Other_Hosts"))

  host_colors <- ifelse(host_groups %in% names(custom_phylum_colors),
                        custom_phylum_colors[host_groups],
                        "#e0e0e0") 
  names(host_colors) <- host_groups

  final_colors <- c(
    host_colors,
    "Virus_Broad"  = "black",
    "Virus_Shared" = "grey60",
    "Other_Hosts"  = "grey85"
  )
  return(final_colors)
}

plot_family_network <- function(output_name) {

  shared_votus <- data %>%
    filter(!is.na(Family), Family != "None") %>%
    group_by(vOTU) %>%
    summarise(n = n_distinct(Family)) %>%
    filter(n > 1) %>% pull(vOTU)

  if(length(shared_votus) == 0) return(NULL)

  edges <- data %>%
    filter(vOTU %in% shared_votus) %>%
    select(vOTU, Host = Family) %>%
    distinct()

  graph <- as_tbl_graph(edges, directed = FALSE) %>%
    activate(nodes) %>%
    left_join(votu_broad_status, by = c("name" = "vOTU")) %>%
    left_join(fam_to_phylum, by = c("name" = "Family")) %>%
    mutate(
      node_type = case_when(
        !is.na(total_species_count) & is_broad == "Broad"  ~ "Virus_Broad",
        !is.na(total_species_count) & is_broad == "Narrow" ~ "Virus_Shared",
        TRUE                                               ~ "Host"
      ),
      color_group = ifelse(node_type == "Host", Phylum, node_type)
    )

  all_groups <- graph %>% as_tibble() %>% pull(color_group) %>% unique()
  my_palette <- get_phylum_palette(all_groups)

  p <- ggraph(graph, layout = 'nicely') +
    geom_edge_link(color = "grey85", width = 0.5, alpha = 0.6) +
    geom_node_point(aes(shape = node_type, color = color_group, size = node_type)) +
    geom_node_text(aes(label = ifelse(node_type != "Virus_Shared", name, "")),
                   repel = TRUE, size = 3, family = "sans", max.overlaps = 50) +
    scale_shape_manual(values = c("Host" = 15, "Virus_Broad" = 16, "Virus_Shared" = 16)) +
    scale_size_manual(values = c("Host" = 5, "Virus_Broad" = 4, "Virus_Shared" = 2)) +
    scale_color_manual(values = my_palette, name = "Phylum / Virus Type") +
    theme_graph(base_family = "sans") +
    labs(title = "Shared vOTU Network (Family Level)",
         subtitle = "Host nodes colored by Phylum (Custom Palette). Black=Broad, Grey=Shared.")

  ggsave(output_name, p, width = 12, height = 10)
}

plot_genus_network <- function(output_name) {

  shared_votus <- data %>%
    filter(!is.na(Genus), Genus != "None") %>%
    group_by(vOTU) %>%
    summarise(n = n_distinct(Genus)) %>%
    filter(n > 1) %>% pull(vOTU)

  edges <- data %>%
    filter(vOTU %in% shared_votus) %>%
    select(vOTU, Host = Genus) %>%
    distinct()

  graph <- as_tbl_graph(edges, directed = FALSE) %>%
    activate(nodes) %>%
    left_join(votu_broad_status, by = c("name" = "vOTU")) %>%
    left_join(genus_to_family, by = c("name" = "Genus")) %>%
    mutate(
      degree = centrality_degree(),
      node_type = case_when(
        !is.na(total_species_count) & is_broad == "Broad"  ~ "Virus_Broad",
        !is.na(total_species_count) & is_broad == "Narrow" ~ "Virus_Shared",
        TRUE                                               ~ "Host"
      ),
      color_group = ifelse(node_type == "Host", Family, node_type),
      show_label = case_when(
        node_type == "Virus_Broad" ~ TRUE,
        node_type == "Host" & rank(-degree) <= 20 ~ TRUE,
        TRUE ~ FALSE
      )
    )

  all_groups <- graph %>% as_tibble() %>% pull(color_group) %>% unique()
  host_groups <- setdiff(all_groups, c("Virus_Broad", "Virus_Shared"))
  n_groups <- length(host_groups)
  host_colors <- pal_igv()(n_groups)
  names(host_colors) <- host_groups
  my_palette <- c(host_colors, "Virus_Broad" = "black", "Virus_Shared" = "grey60")

  p <- ggraph(graph, layout = 'nicely') +
    geom_edge_link(color = "grey85", width = 0.4, alpha = 0.5) +
    geom_node_point(aes(shape = node_type, color = color_group, size = node_type)) +
    geom_node_text(aes(label = ifelse(show_label, name, "")),
                   repel = TRUE, size = 3, family = "sans", fontface = "italic", max.overlaps = 50) +
    scale_shape_manual(values = c("Host" = 15, "Virus_Broad" = 16, "Virus_Shared" = 16)) +
    scale_size_manual(values = c("Host" = 4, "Virus_Broad" = 3.5, "Virus_Shared" = 2)) +
    scale_color_manual(values = my_palette, name = "Family / Virus Type") +
    theme_graph(base_family = "sans") +
    labs(title = "Shared vOTU Network (Genus Level)")

  ggsave(output_name, p, width = 14, height = 12)
}

plot_species_network <- function(color_by = "Phylum", output_name) {
  message(paste0("\n>>> 正在绘制: Species Level (Color by ", color_by, ")..."))

  shared_votus_list <- data %>%
    filter(!is.na(Species), Species != "None") %>%
    group_by(vOTU) %>%
    summarise(n_targets = n_distinct(Species)) %>%
    filter(n_targets > 1) %>% pull(vOTU)

  if(length(shared_votus_list) == 0) return(NULL)

  edges <- data %>%
    filter(vOTU %in% shared_votus_list) %>%
    select(vOTU, Species) %>%
    distinct()

  graph_ready <- as_tbl_graph(edges, directed = FALSE) %>%
    activate(nodes) %>%
    left_join(votu_broad_status, by = c("name" = "vOTU")) %>%
    left_join(host_taxonomy, by = c("name" = "Species")) %>%
    mutate(
      degree = centrality_degree(),
      node_class = case_when(
        !is.na(total_species_count) & is_broad == "Broad"  ~ "Virus_Broad",
        !is.na(total_species_count) & is_broad == "Narrow" ~ "Virus_Shared",
        TRUE                                               ~ "Host"
      ),
      show_label = case_when(
        node_class == "Virus_Broad" ~ TRUE,
        node_class == "Host" & rank(-degree) <= 20 ~ TRUE,
        TRUE ~ FALSE
      )
    )

  if (color_by == "Phylum") {

    graph_ready <- graph_ready %>%
      mutate(color_group = ifelse(node_class == "Host", Phylum, node_class))
    
    all_groups <- graph_ready %>% as_tibble() %>% pull(color_group) %>% unique()
    final_colors <- get_phylum_palette(all_groups)

  } else {

    graph_ready <- graph_ready %>%
      mutate(
        family_lumped = fct_lump_n(Family, n = 15, other_level = "Other_Hosts"),
        color_group = ifelse(node_class == "Host", as.character(family_lumped), node_class)
      )
    all_groups <- graph_ready %>% as_tibble() %>% pull(color_group) %>% unique()
    host_groups <- setdiff(all_groups, c("Virus_Broad", "Virus_Shared", "Other_Hosts"))
    base_colors <- scales::hue_pal()(length(host_groups))
    names(base_colors) <- host_groups
    final_colors <- c(base_colors, "Virus_Broad"="black", "Virus_Shared"="grey60", "Other_Hosts"="grey85")
  }

  p <- ggraph(graph_ready, layout = 'stress') +
    geom_edge_link(color = "grey85", width = 0.3, alpha = 0.5) +
    geom_node_point(aes(shape = node_class, color = color_group, size = node_class)) +
    geom_node_text(aes(label = ifelse(show_label, name, "")),
                   repel = TRUE, size = 3, family = "sans", fontface = "italic", max.overlaps = 50) +
    scale_shape_manual(values = c("Host" = 15, "Virus_Broad" = 16, "Virus_Shared" = 16)) +
    scale_size_manual(values = c("Host" = 3, "Virus_Broad" = 3, "Virus_Shared" = 1.5)) +
    scale_color_manual(values = final_colors, name = "Classification") +
    theme_graph(base_family = "sans") +
    labs(title = paste0("Species Network - Color by ", color_by))

  ggsave(output_name, p, width = 14, height = 12)
}

plot_family_network("Network_Family_CustomPhylum.pdf")
plot_genus_network("Network_Genus_FamilyColor.pdf")
plot_species_network(color_by = "Phylum", "Network_Species_CustomPhylum.pdf")
plot_species_network(color_by = "Family", "Network_Species_FamilyTop15.pdf")
