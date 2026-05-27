
library(data.table)
library(tidyverse)
library(ggplot2)
library(scales) 

df <- fread("active_prophage_stats.tsv")

uhgv_cols <- c("v_Family", "v_SubFamily", "v_Genus", "v_SubGenus", "v_MemberID")
ictv_cols <- c("ictv_Realm", "ictv_Kingdom", "ictv_Phylum", "ictv_Class", "ictv_Order", "ictv_Family", "ictv_Genus")

get_lowest_rank <- function(data, cols, prefix) {
  res <- rep("Unclassified", nrow(data))

  for (col in cols) {
    if (col %in% names(data)) {
      val <- data[[col]]
      valid_idx <- which(!is.na(val) & trimws(val) != "" & tolower(trimws(val)) != "unclassified")
      clean_name <- gsub(prefix, "", col)
      res[valid_idx] <- clean_name
    }
  }
  return(res)
}

df$UHGV_Depth <- get_lowest_rank(df, uhgv_cols, "v_")
df$ICTV_Depth <- get_lowest_rank(df, ictv_cols, "ictv_")

uhgv_summary <- df %>% 
  count(UHGV_Depth) %>% 
  mutate(Database = "UHGV Lowest Rank",
         prop = n / sum(n)) %>%
  rename(Rank = UHGV_Depth)

ictv_summary <- df %>% 
  count(ICTV_Depth) %>% 
  mutate(Database = "ICTV Lowest Rank",
         prop = n / sum(n)) %>%
  rename(Rank = ICTV_Depth)

combined_summary <- bind_rows(uhgv_summary, ictv_summary)

rank_levels <- c("Realm", "Kingdom", "Phylum", "Class", "Order", 
                 "Family", "SubFamily", "Genus", "SubGenus", "MemberID", "Unclassified")
combined_summary$Rank <- factor(combined_summary$Rank, levels = intersect(rank_levels, unique(combined_summary$Rank)))

unique_ranks <- levels(combined_summary$Rank)
my_colors <- scales::hue_pal()(length(unique_ranks))
names(my_colors) <- unique_ranks
if ("Unclassified" %in% names(my_colors)) {
  my_colors["Unclassified"] <- "#E0E0E0" 
}

p <- ggplot(combined_summary, aes(x = "", y = prop, fill = Rank)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~ Database) +
  geom_text(aes(label = ifelse(prop > 0.02, scales::percent(prop, accuracy = 0.1), "")),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = my_colors) +
  theme_void() +
  theme(
    strip.text = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right"
  )


ggsave("taxonomy_annotation_depth.pdf", p, width = 10, height = 5)

