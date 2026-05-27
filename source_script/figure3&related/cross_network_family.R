
suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
  library(pheatmap)
  library(ggsci)
  library(reshape2)
  library(igraph)
  library(ggraph)
  library(tidygraph)
})


# setwd("~/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats")


file_path <- "phage_full_stats_separated.tsv"

data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df_clean <- data %>%
  filter(!is.na(vOTU) & vOTU != "") %>%
  filter(!is.na(Family) & Family != "") %>%  
  select(vOTU, Family, Phylum) %>%           
  distinct() 

votu_counts <- df_clean %>%
  count(vOTU) %>%
  filter(n > 1) 
cross_family_df <- df_clean %>%
  filter(vOTU %in% votu_counts$vOTU)


if(nrow(cross_family_df) == 0) stop("no cross family vOTU")

host_phage_mat <- table(cross_family_df$Family, cross_family_df$vOTU)
host_co_mat <- as.matrix(host_phage_mat %*% t(host_phage_mat))
diag(host_co_mat) <- 0

edge_list <- as.data.frame(as.table(host_co_mat)) %>%
  filter(Freq > 0) %>%
  rename(From = Var1, To = Var2, Weight = Freq) %>%
  filter(as.character(From) < as.character(To)) %>%
  mutate(From = as.character(From),
         To = as.character(To)) 
write.csv(edge_list, "Family_Network_Edges_Summary.csv", row.names = FALSE) 


edge_list_filtered <- edge_list %>% filter(Weight >= 2)
node_info <- df_clean %>%
  select(Family, Phylum) %>%                
  distinct(Family, .keep_all = TRUE) %>%    
  mutate(Family = as.character(Family))     
graph <- as_tbl_graph(edge_list, directed = FALSE) %>%
  activate(nodes) %>%
  left_join(node_info, by = c("name" = "Family")) %>%  
  mutate(degree = centrality_degree(weights = Weight)) %>%
  filter(!node_is_isolated())

pdf("1_Family_Network_Topology.pdf", width = 12, height = 10) 
set.seed(123) 

p <- ggraph(graph, layout = "graphopt", charge = 0.05, spring.length = 100) +
  geom_edge_link(aes(width = Weight, alpha = Weight), color = "black") +
  scale_edge_width(range = c(0.5, 5), name = "Shared vOTUs") +
  scale_edge_alpha(range = c(0.2, 0.95), guide = "none") +
  geom_node_point(aes(size = degree, fill = Phylum), shape = 21, color = "black", stroke = 0.5) +
  scale_size(range = c(2, 20), name = "Connectivity") +
  scale_fill_npg(name = "Host Phylum") + 

  geom_node_text(aes(label = ifelse(degree > quantile(degree, 0.8, na.rm=TRUE), name, "")),
                 repel = TRUE, size = 3.5, fontface = "bold", bg.color = "white", bg.r = 0.15) +

  theme_graph(base_family = "sans") +
  labs(title = "Host Family Connectivity Network",  
       subtitle = paste0("Nodes: ", length(unique(c(edge_list$From, edge_list$To))),
                         " Families | Edges: ", nrow(edge_list), " Connections"), 
       caption = "Edge width correlates with number of shared vOTUs")

print(p)
dev.off()

