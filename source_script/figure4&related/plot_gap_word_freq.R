library(tidyverse)
library(patchwork) 

df <- read.csv("packaging_module_structure_vFAM7.csv")

clean_tags <- function(column_data) {
  exclude <- c("hypothetical protein", "unknown function", "None", 
               "Missing_A", "Missing_C", "Missing_TerS", "Missing_TerL", "Unknown", "")
  
  column_data %>%
    as.character() %>%
    str_split(";") %>%          
    unlist() %>%                 
    str_trim() %>%               
    .[!(. %in% exclude)] %>%     
    table() %>%                  
    as.data.frame() %>%
    rename(Feature = 1, Count = 2) %>%
    arrange(desc(Count)) %>%
    head(15)                     
}

ab_annot <- clean_tags(df$Gap_HNH_TerS_Annot)
ab_cat   <- clean_tags(df$Gap_HNH_TerS_Category)
bc_annot <- clean_tags(df$Gap_TerS_TerL_Annot)
bc_cat   <- clean_tags(df$Gap_TerS_TerL_Category)
plot_freq <- function(data, title, fill_color) {
  ggplot(data, aes(x = reorder(Feature, Count), y = Count)) +
    geom_bar(stat = "identity", fill = fill_color, alpha = 0.8) +
    coord_flip() + # 横向条形图，方便阅读长基因名
    theme_minimal() +
    labs(title = title, x = "", y = "Occurrences") +
    theme(plot.title = element_text(size = 10, face = "bold"))
}


p_ab <- plot_freq(ab_annot, "Top Annotations (AB Gap)", "steelblue") +
        plot_freq(ab_cat, "Top Categories (AB Gap)", "darkblue") +
        plot_annotation(title = "Functional Content in HNH - TerS Intergenic Region")


p_bc <- plot_freq(bc_annot, "Top Annotations (BC Gap)", "darkorange") +
        plot_freq(bc_cat, "Top Categories (BC Gap)", "darkred") +
        plot_annotation(title = "Functional Content in TerS - TerL Intergenic Region")


ggsave("Freq_AB_Content.pdf", p_ab, width = 12, height = 6)
ggsave("Freq_BC_Content.pdf", p_bc, width = 12, height = 6)

print("Visualization complete: Freq_AB_Content.pdf and Freq_BC_Content.pdf")

