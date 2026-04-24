
library(tidyverse)
library(scales)

BAR_COLOR_LOW  <- "grey85"   
BAR_COLOR_HIGH <- "#447095"  
BAR_BORDER_COLOR <- "grey30" 
LABEL_TEXT_SIZE  <- 5      
LABEL_TEXT_COLOR <- "black"  
LABEL_FONT_FACE  <- "bold"   
LABEL_VJUST      <- -0.5     
AXIS_TITLE_SIZE  <- 16       
AXIS_TITLE_COLOR <- "black"  
AXIS_TEXT_SIZE   <- 14       
AXIS_TEXT_COLOR  <- "black"  
PLOT_TITLE_SIZE  <- 16       
PLOT_TITLE_COLOR <- "black"  
SUBTITLE_SIZE    <- 12       
LEGEND_TITLE_SIZE <- 12      
LEGEND_TEXT_SIZE  <- 11      
LEGEND_POSITION   <- "right" 
PLOT_WIDTH  <- 12
PLOT_HEIGHT <- 5

file_name <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/sig_act_prevalence/clustering_results/vfc_full_flow_breakdown.tsv"

df <- read_tsv(file_name, show_col_types = FALSE)

plot_data <- df %>%
  mutate(
    VFC_Label = paste0("VFC_", str_remove(as.character(VFC_Source), "\\.0$")),
    User_Significant_Count = as.numeric(User_Significant_Count)
  ) %>%
  group_by(VFC_Source, VFC_Label) %>%
  summarize(
    Total_Significant_Count = sum(User_Significant_Count, na.rm = TRUE),
    .groups = "drop" 
  ) %>%
  arrange(VFC_Source) %>%
  mutate(VFC_Label = fct_inorder(VFC_Label))

p_bar <- ggplot(plot_data, aes(x = VFC_Label, y = Total_Significant_Count, fill = Total_Significant_Count)) +
  geom_col(color = BAR_BORDER_COLOR, width = 0.7, alpha = 0.9) +
  geom_text(aes(label = ifelse(Total_Significant_Count > 0, Total_Significant_Count, "")),
            color = LABEL_TEXT_COLOR,
            size = LABEL_TEXT_SIZE,
            fontface = LABEL_FONT_FACE,
            vjust = LABEL_VJUST) +
  scale_fill_gradient(low = BAR_COLOR_LOW, high = BAR_COLOR_HIGH,
                      name = "Significant\nvOTU Count") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  
  labs(
    title = "Total Significant User Prophage Activity by VFC Source",
    subtitle = paste0("Bars represent the sum of significant vOTUs originating from each Reference Cluster (VFC)"),
    x = "Reference Cluster (VFC)",
    y = "Total Significant vOTU Count"
  ) +
  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = AXIS_TITLE_SIZE, color = AXIS_TITLE_COLOR, face = "bold"),
    axis.title.y = element_text(size = AXIS_TITLE_SIZE, color = AXIS_TITLE_COLOR, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = AXIS_TEXT_SIZE, color = AXIS_TEXT_COLOR),
    axis.text.y = element_text(size = AXIS_TEXT_SIZE, color = AXIS_TEXT_COLOR),
    
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(linetype = "dashed", color = "grey90"),
    panel.grid.minor = element_blank(),
    
    legend.position = LEGEND_POSITION,
    legend.title = element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
    legend.text = element_text(size = LEGEND_TEXT_SIZE),
    
    plot.title = element_text(size = PLOT_TITLE_SIZE, color = PLOT_TITLE_COLOR, face = "bold"),
    plot.subtitle = element_text(size = SUBTITLE_SIZE)
  )


ggsave("VFC_Bar_Plot.pdf", p_bar, width = PLOT_WIDTH, height = PLOT_HEIGHT)



