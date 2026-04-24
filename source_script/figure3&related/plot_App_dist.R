library(data.table)
library(tidyverse)
library(ggplot2)

process_dist_file <- function(file_path, db_name) {

  dt <- fread(file_path, select = c(1, 3), header = FALSE, col.names = c("query_raw", "dist"))
  dt[, query_name := basename(query_raw)]

  dt_min <- dt[, .(min_dist = min(dist, na.rm = TRUE)), by = query_name]
  setnames(dt_min, "min_dist", db_name)
  
  return(dt_min)
}

df_inphared   <- process_dist_file("active_votu-inphared_dist.txt", "inphared")
df_prophagedb <- process_dist_file("active_votu-prophagedb_dist.txt", "prophagedb")
df_uhgv       <- process_dist_file("active_votu-uhgv_dist.txt", "uhgv")

merged_df <- df_inphared %>%
  full_join(df_prophagedb, by = "query_name") %>%
  full_join(df_uhgv, by = "query_name")

print(head(merged_df))
fwrite(merged_df, "min_genetic_distances_merged.tsv", sep = "\t", na = "NA")

long_df <- merged_df %>%
  pivot_longer(cols = c("inphared", "prophagedb", "uhgv"),
               names_to = "Database",
               values_to = "Genetic_Distance") %>%
  filter(!is.na(Genetic_Distance)) %>% # 去除空值以便绘图

  mutate(Database_fct = factor(Database, levels = c("inphared", "prophagedb", "uhgv")),
         x_num = as.numeric(Database_fct))

p <- ggplot(long_df, aes(y = Genetic_Distance, fill = Database_fct)) +
  geom_jitter(aes(x = x_num + 0.15), 
              width = 0.1, size = 0.8, alpha = 0.3, color = "grey40") +
  geom_boxplot(aes(x = x_num - 0.15, group = x_num), 
               width = 0.2, alpha = 0.8, outlier.shape = NA, color = "black") +
  scale_x_continuous(breaks = 1:3, labels = levels(long_df$Database_fct)) +
  theme_bw() +
  labs(title = "Minimum Genetic Distance Distribution",
       subtitle = "Split Boxplot and Scatter comparison across databases",
       x = "Database",
       y = "Min Genetic Distance") +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14),
        legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

ggsave("genetic_distance_plot_split.pdf", p, width = 8, height = 6)

