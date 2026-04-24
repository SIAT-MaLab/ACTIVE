
suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

ACTIVE_THRESHOLD <- 0.7
MIN_GENOME_COUNT <- 5   

cross_colors <- c(
  "Non-crossing"  = "grey80", "Cross-Species" = "#A8DADC", "Cross-Genus"   = "#457B9D",
  "Cross-Family"  = "#E63946", "Cross-Order"   = "#D62828", "Cross-Class"   = "#9D0208", "Cross-Phylum"  = "#000000"
)
level_order <- c("Non-crossing", "Cross-Species", "Cross-Genus", "Cross-Family", "Cross-Order", "Cross-Class", "Cross-Phylum")
group_colors <- c("Inactive vOTU" = "#BDC3C7", "Other Active vOTU" = "#3498DB", "SAvOTU" = "#E74C3C")
group_levels <- c("Inactive vOTU", "Other Active vOTU", "SAvOTU")

PHYLUM_COLORS <- c(
  "Bacillota_A" = "#dce8ef", "Bacillota_C" = "#add2e5", "Bacillota_I" = "#34666b", "Bacillota" = "#7abbce",
  "Bacillota_B" = "#769499", "Actinomycetota" = "#f5cfa6", "Bacteroidota" = "#efe3ef", "Verrucomicrobiota" = "#d4d68a",
  "Pseudomonadota" = "#e68f9f", "Desulfobacterota" = "#d5d5d6", "Fusobacteriota" = "#b7a3c9"
)
STROKE_WIDTH <- 0.8
STROKE_COLOR_DEFAULT <- "darkgray"
STROKE_COLOR_CROSS   <- "black"

raw_meta <- read_tsv("merged_phage_stats_taxonomy.tsv", show_col_types = FALSE)

df_clean <- raw_meta %>%
  mutate(Activity_Score = suppressWarnings(as.numeric(as.character(Activity_Score)))) %>%
  mutate(Activity_Score = replace_na(Activity_Score, 0)) %>%
  mutate(Is_Active = ifelse(Activity_Score >= ACTIVE_THRESHOLD, "Yes", "No")) %>%
  mutate(across(c(Species, Genus, Family, Order, Class, Phylum), 
                ~ na_if(str_trim(.), "") %>% na_if("s__") %>% na_if("g__") %>% na_if("f__") %>%
                  na_if("Unclassified") %>% na_if("unknown")))

get_mode <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if(length(x) == 0) return(NA)
  names(table(x))[which.max(table(x))]
}


run_fisher <- function(data, level_name) {
  data_lvl <- data %>% filter(!is.na(!!sym(level_name)) & !!sym(level_name) != "")
  tot_active <- sum(data_lvl$Is_Active == "Yes")
  tot_non_active <- sum(data_lvl$Is_Active == "No")
  
  data_lvl %>%
    group_by(!!sym(level_name)) %>%
    summarise(
      Count_Active = sum(Is_Active == "Yes"),
      Count_Total = n(),
      Prevalence_Pct = (Count_Active / Count_Total) * 100,
      Mean_Intensity_When_Active = mean(Activity_Score[Activity_Score >= ACTIVE_THRESHOLD]),
      Max_Intensity = max(Activity_Score),
      Phylum = get_mode(Phylum),
      Genome_Count = n_distinct(Genome),
      .groups = "drop"
    ) %>%
    mutate(Mean_Intensity_When_Active = replace_na(Mean_Intensity_When_Active, 0)) %>%
    filter(Count_Active > 0, Count_Total >= MIN_GENOME_COUNT) %>%
    rowwise() %>%
    mutate(
      p_value = fisher.test(matrix(c(Count_Active, tot_active - Count_Active, 
                                     Count_Total - Count_Active, tot_non_active - (Count_Total - Count_Active)), 
                                   nrow=2), alternative = "greater")$p.value
    ) %>%
    ungroup() %>%
    mutate(
      p.adj = p.adjust(p_value, method = "fdr"),
      Significance = ifelse(p.adj < 0.05, "Significant", "Not Significant")
    ) %>% arrange(desc(Prevalence_Pct))
}

fisher_votu <- run_fisher(df_clean, "vOTU")
fisher_vfam <- run_fisher(df_clean, "vFAM")

write_tsv(fisher_votu, "stats_fisher_vOTU.tsv")
write_tsv(fisher_vfam, "stats_fisher_vFAM.tsv")


votu_cross_stats <- df_clean %>%
  filter(!is.na(vOTU) & vOTU != "") %>%
  group_by(vOTU) %>%
  summarise(
    Member_Count = n(),
    n_Species = n_distinct(Species, na.rm = TRUE),
    n_Genus   = n_distinct(Genus, na.rm = TRUE),
    n_Family  = n_distinct(Family, na.rm = TRUE),
    n_Order   = n_distinct(Order, na.rm = TRUE),
    n_Class   = n_distinct(Class, na.rm = TRUE),
    n_Phylum  = n_distinct(Phylum, na.rm = TRUE),
    Max_Activity = max(Activity_Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Cross_Category = case_when(
      n_Phylum > 1 ~ "Cross-Phylum", n_Class > 1 ~ "Cross-Class", n_Order > 1 ~ "Cross-Order",
      n_Family > 1 ~ "Cross-Family", n_Genus > 1 ~ "Cross-Genus", n_Species > 1 ~ "Cross-Species",
      TRUE         ~ "Non-crossing"
    ),
    Cross_Category = factor(Cross_Category, levels = level_order),
    Is_Generalist = if_else(n_Species > 1, "Generalist", "Specialist")
  )

df_unified <- votu_cross_stats %>%
  left_join(fisher_votu %>% select(vOTU, Significance), by = "vOTU") %>%
  mutate(
    Significance = replace_na(Significance, "Not Significant"), 
    Analysis_Group = case_when(
      Max_Activity >= 0.7 & Significance == "Significant" ~ "SAvOTU",
      Max_Activity >= 0.7 ~ "Other Active vOTU",
      TRUE ~ "Inactive vOTU" 
    ),
    Analysis_Group = factor(Analysis_Group, levels = group_levels)
  )

df_full_detailed <- df_unified %>%
  select(vOTU, Member_Count, Cross_Category, Is_Generalist, Analysis_Group, Max_Activity, n_Phylum:n_Species) %>%
  arrange(desc(Cross_Category), desc(Member_Count))

write_csv(df_full_detailed, "Detailed_Cross_Host_vOTUs.csv")


draw_bubble <- function(fisher_data, level_name) {
  if(level_name == "vOTU") {
    plot_df <- fisher_data %>%
      left_join(df_full_detailed %>% select(vOTU, Cross_Category), by = "vOTU") %>%
      mutate(
        Is_Cross_Host = !is.na(Cross_Category) & Cross_Category != "Non-crossing",
        Cross_Label = case_when(
          Significance == "Not Significant" ~ "No Cross / Not Sig",
          Is_Cross_Host ~ "Cross-Host",
          TRUE ~ "No Cross / Not Sig"
        )
      )
  } else {
    plot_df <- fisher_data %>% mutate(Cross_Label = "No Cross / Not Sig")
  }
  
  plot_df <- plot_df %>%
    mutate(
      Color_Group = ifelse(Significance == "Significant", Phylum, "Not Significant"),
      Plot_Order = ifelse(Significance == "Significant", 1, 0),
      log_FDR = -log10(ifelse(p.adj == 0, .Machine$double.xmin, p.adj))
    ) %>% arrange(Plot_Order)
    
  actual_phyla <- unique(plot_df$Phylum[plot_df$Significance == "Significant"])
  actual_phyla <- actual_phyla[!is.na(actual_phyla)]
  current_palette <- setNames(rep("#bdbdbd", length(actual_phyla)), actual_phyla)
  for(p in actual_phyla) { if(p %in% names(PHYLUM_COLORS)) current_palette[p] <- PHYLUM_COLORS[p] }
  final_colors <- c(current_palette, "Not Significant" = "grey95")
  
  top_hits <- plot_df %>% arrange(p.adj) %>% head(10)
  sig_threshold <- -log10(0.05) 
  
  p <- ggplot(plot_df, aes(x = log_FDR, y = Prevalence_Pct)) +
    geom_vline(xintercept = sig_threshold, linetype = "dashed", color = "red", alpha = 0.6, linewidth = 0.8) +
    geom_point(aes(size = Genome_Count, fill = Color_Group, color = Cross_Label), stroke = STROKE_WIDTH, shape = 21, alpha = 1.0) +
    scale_fill_manual(values = final_colors, name = "Host Phylum") +
    scale_color_manual(name = "Cross-Host Status", values = c("Cross-Host" = STROKE_COLOR_CROSS, "No Cross / Not Sig" = STROKE_COLOR_DEFAULT)) +
    scale_size_continuous(range = c(2, 12), name = "Genome Count") +
    labs(title = paste0("Activity Landscape: ", level_name),
         subtitle = paste0("Fisher Test Significant (FDR < 0.05): n=", sum(plot_df$Significance == "Significant")),
         x = "-log10(FDR) [Statistical Significance]", y = "% of Genome Active") +
    theme_minimal() +
    theme(aspect.ratio = 0.75, plot.title = element_text(size=16, face="bold"), panel.border = element_rect(colour="black", fill=NA)) +
    guides(fill = guide_legend(override.aes = list(shape=21, size=5, color=STROKE_COLOR_DEFAULT)))
  
  if(nrow(top_hits) > 0) {
    p <- p + geom_text_repel(data = top_hits, aes(label = !!sym(level_name)), size=4.5, fontface="bold",
                             min.segment.length=0, segment.color="black", box.padding=0.8, max.overlaps=Inf)
  }
  return(p)
}

bubble_plot <- draw_bubble(fisher_votu, "vOTU") / draw_bubble(fisher_vfam, "vFAM") + plot_annotation(tag_levels = 'A')
ggsave("Combined_Activity_Landscape_vOTU_vFAM.pdf", bubble_plot, width = 12, height = 18, bg = "white")

df_A <- df_unified %>% filter(Analysis_Group == "SAvOTU") %>% count(Cross_Category) %>% mutate(Prop = n/sum(n))
write_csv(df_A, "Table_A_SAvOTU_Cross_Host_Categories_Summary.csv")

pA <- ggplot(df_A, aes(x = Cross_Category, y = n, fill = Cross_Category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, Prop*100)), vjust = -0.2, size = 3.5) +
  scale_fill_manual(values = cross_colors) + scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "A. Cross-Host Distribution of SAvOTUs", x = "", y = "Number of SAvOTUs") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1, face="bold"))

df_B <- df_unified %>% group_by(Analysis_Group) %>%
  summarise(Total_vOTU = n(), Generalists = sum(Is_Generalist == "Generalist"), Prop_Generalist = Generalists/Total_vOTU, .groups="drop")

calc_fisher <- function(g1, g2) { fisher.test(table((df_unified %>% filter(Analysis_Group %in% c(g1, g2)))$Analysis_Group, 
                                                    (df_unified %>% filter(Analysis_Group %in% c(g1, g2)))$Is_Generalist))$p.value }
pval_S_O <- calc_fisher("SAvOTU", "Other Active vOTU")
pval_S_I <- calc_fisher("SAvOTU", "Inactive vOTU")
write_csv(df_B, "Table_B_Generalist_Proportion.csv")

format_pval <- function(p) { if(p < 0.0001) "****" else if(p < 0.001) "***" else if(p < 0.01) "**" else if(p < 0.05) "*" else "ns" }
ymax <- max(df_B$Prop_Generalist)

pB <- ggplot(df_B, aes(x = Analysis_Group, y = Prop_Generalist, fill = Analysis_Group)) +
  geom_bar(stat="identity", width=0.6, color="black") +
  geom_text(aes(label = percent(Prop_Generalist, accuracy=0.1)), vjust=-0.5, size=4) +
  geom_segment(aes(x=2, xend=3, y=ymax*1.15, yend=ymax*1.15)) +
  geom_text(aes(x=2.5, y=ymax*1.18), label=sprintf("%s\n(p=%.2e)", format_pval(pval_S_O), pval_S_O), size=3) +
  geom_segment(aes(x=1, xend=3, y=ymax*1.35, yend=ymax*1.35)) +
  geom_text(aes(x=2, y=ymax*1.38), label=sprintf("%s\n(p=%.2e)", format_pval(pval_S_I), pval_S_I), size=3) +
  scale_fill_manual(values = group_colors) + scale_y_continuous(labels = percent, limits = c(0, ymax*1.5)) +
  labs(title = "B. Proportion of Generalists", x="", y="Generalist Proportion (%)") +
  theme_classic() + theme(legend.position="none", axis.text.x = element_text(angle=15, hjust=1, face="bold"))

comprehensive_plot <- pA / pB + plot_annotation(title = "SAvOTU vs Other Active vs Inactive vOTUs")

ggsave("Fig_Comprehensive_2Groups.pdf", comprehensive_plot, width = 8, height = 9)

