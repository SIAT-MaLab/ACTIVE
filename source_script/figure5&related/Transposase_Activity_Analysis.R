

library(tidyverse)
library(patchwork)
library(ggpubr)
library(scales)

metadata_file  <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_checkv_quality/hqphage_metadata.tsv"
lifestyle_file <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_annot/pfamscan/Phage_Lifestyle.tsv"

cat("Loading data...\n")
meta_df <- read_tsv(metadata_file, show_col_types = FALSE)
life_df <- read_tsv(lifestyle_file, show_col_types = FALSE)


format_pval <- function(p) {
  if(p < 0.001) "P < 0.001" else paste0("P = ", signif(p, 3))
}

clean_df <- meta_df %>%
  inner_join(life_df, by = c("Phage_Name" = "Contig")) %>%
  filter(!is.na(ictv_Class) & grepl("Caudoviricetes", ictv_Class)) %>%
  mutate(
    is_active = !is.na(Activity_Score) & Activity_Score >= 0.7,
    Group_Tnp = factor(ifelse(Encodes_Transposase == "Yes", "Transposase (+)", "Transposase (-)"), 
                       levels = c("Transposase (+)", "Transposase (-)")),
    Group_Mu  = factor(ifelse(Analyzed_Mu_Lifestyle == "Mu-like Lifestyle", "Mu-like", "Other"), 
                       levels = c("Mu-like", "Other"))
  )

write_tsv(clean_df, "clean_df_output.tsv")
cat("Clean data exported to clean_df_output.tsv\n")
# ==============================================================================
# TASK 2: Encodes_Transposase 
# ==============================================================================
cat("\n--- Running Task 2: Transposase Encoding Analysis (Overall Only) ---\n")

stats_t2 <- clean_df %>% 
  group_by(Group_Tnp) %>% 
  summarise(Total = n(), Active = sum(is_active), Rate = Active/Total, .groups = "drop")
f_test_t2 <- fisher.test(table(clean_df$Group_Tnp, clean_df$is_active))

p_t2_main <- ggplot(stats_t2, aes(x = Group_Tnp, y = Rate, fill = Group_Tnp)) +
  geom_col(width = 0.5, color = "black") +
  geom_text(aes(label = paste0(round(Rate*100, 1), "%\n(n=", Total, ")")), vjust = -0.3) +
  scale_y_continuous(labels = percent, limits = c(0, max(stats_t2$Rate) * 1.25)) +
  scale_fill_manual(values = setNames(c("#E74C3C", "#95A5A6"), c("Transposase (+)", "Transposase (-)"))) +
  labs(title = "Task 2: Overall Activity by Transposase", 
       subtitle = paste0("Fisher Test: ", format_pval(f_test_t2$p.value)), 
       x = "", y = "Active Rate") +
  theme_classic() + theme(legend.position = "none")

ggsave("Task2_Transposase_Overall_Activity.pdf", p_t2_main, width = 6, height = 5)

# ==============================================================================
# TASK 3: Mu-like Lifestyle 
# ==============================================================================
cat("\n--- Running Task 3: Mu-like Lifestyle Analysis ---\n")


stats_t3_main <- clean_df %>% group_by(Group_Mu) %>% 
  summarise(Total = n(), Active = sum(is_active), Rate = Active/Total, .groups = "drop")
f_test_t3 <- fisher.test(table(clean_df$Group_Mu, clean_df$is_active))

p_t3_main <- ggplot(stats_t3_main, aes(x = Group_Mu, y = Rate, fill = Group_Mu)) +
  geom_col(width = 0.6, color = "black") +
  geom_text(aes(label = paste0(round(Rate*100, 1), "%\n(n=", Total, ")")), vjust = -0.3) +
  scale_y_continuous(labels = percent, limits = c(0, max(stats_t3_main$Rate) * 1.25)) +
  scale_fill_manual(values = setNames(c("#9B59B6", "#95A5A6"), c("Mu-like", "Other"))) +
  labs(title = "A. Overall Mu-like Activity", subtitle = paste0("Fisher: ", format_pval(f_test_t3$p.value)), x = "", y = "Active Rate") +
  theme_classic() + theme(legend.position = "none")


cat("Generating vFAM summary table...\n")

global_vfam_counts <- clean_df %>% 
  filter(!is.na(vFAM) & vFAM != "NA" & vFAM != "") %>%
  count(vFAM, name = "Global_Total")

mu_phages <- clean_df %>% filter(Group_Mu == "Mu-like")

vfam_summary <- mu_phages %>%
  filter(!is.na(vFAM) & vFAM != "NA" & vFAM != "") %>%
  group_by(vFAM) %>%
  summarise(
    Total_Members = n(), 
    Active_Count = sum(is_active),
    Active_Rate = Active_Count / Total_Members,
    .groups = "drop"
  ) %>%
  left_join(global_vfam_counts, by = "vFAM") %>%
  mutate(Mu_Proportion = Total_Members / Global_Total) %>% # 计算 Mu-like 占比 (% in vFAM)
  arrange(desc(Total_Members))

write_tsv(vfam_summary, "Task3_Mu_like_vFAM_Summary.tsv")

vfam_plot_df <- vfam_summary %>% 
  filter(Total_Members > 10) %>%
  mutate(vFAM = factor(vFAM, levels = rev(vFAM))) 

p_b1 <- ggplot(vfam_plot_df, aes(x = vFAM, y = Total_Members)) + 
  geom_col(fill = "#3498DB", width = 0.7) +
  geom_text(aes(label = Total_Members), hjust = -0.2, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + coord_flip() +
  labs(title = "B1. vFAM Count (n > 10)", x = "", y = "Count of Mu-like") + theme_bw()

p_b2 <- ggplot(vfam_plot_df, aes(x = vFAM, y = Active_Rate)) + 
  geom_col(fill = "#5DADE2", color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(Active_Rate*100,0),"%")), hjust = -0.2, size = 3) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.3))) + coord_flip() +
  labs(title = "B2. Active Rate", x = "", y = "Rate") +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_b3 <- ggplot(vfam_plot_df, aes(x = vFAM, y = Mu_Proportion)) + 
  geom_col(fill = "#8E44AD", alpha=0.7, width = 0.7) +
  geom_text(aes(label = paste0(round(Mu_Proportion*100,0),"%")), hjust = -0.2, size = 3) +
  scale_y_continuous(labels = percent, limits = c(0, 1.15), expand = expansion(mult = c(0, 0))) + coord_flip() +
  labs(title = "B3. % in vFAM", subtitle = "Prop. of Mu-like", x = "", y = "Proportion") +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

layout_b <- p_b1 + p_b2 + p_b3 + plot_layout(widths = c(1.2, 1, 1))


cat("Generating strict host family plots...\n")

strict_vfams <- vfam_summary %>%
  filter(Total_Members > 10 & Mu_Proportion >= 0.8) %>%
  pull(vFAM)

host_plot_df <- mu_phages %>%
  filter(vFAM %in% strict_vfams) %>%
  filter(!is.na(Family) & Family != "Unknown" & Family != "NA" & Family != "") %>%
  group_by(Family) %>%
  summarise(Total_Members = n(), Active_Rate = sum(is_active)/n(), .groups = "drop") %>%
  arrange(desc(Total_Members)) %>%
  mutate(Family = factor(Family, levels = rev(Family)))

if(nrow(host_plot_df) > 0) {
  p_c1 <- ggplot(host_plot_df, aes(x = Family, y = Total_Members)) + 
    geom_col(fill = "#2ECC71", width = 0.7) +
    geom_text(aes(label = Total_Members), hjust = -0.2, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + coord_flip() +
    labs(title = "C1. Host Count (Strict vFAMs)", x = "", y = "Count") + theme_bw()
  
  p_c2 <- ggplot(host_plot_df, aes(x = Family, y = Active_Rate)) + 
    geom_col(fill = "#58D68D", color = "black", width = 0.7) +
    geom_text(aes(label = paste0(round(Active_Rate*100,0),"%")), hjust = -0.2, size = 3) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.3))) + coord_flip() +
    labs(title = "C2. Active Rate", x = "", y = "Rate") +
    theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  layout_c <- p_c1 + p_c2 + plot_layout(widths = c(1, 1))
} else {
  layout_c <- ggplot() + theme_void() + 
    labs(title = "No data passed the strict C filter\n(Member>10 & Mu_Prop>=80%)")
}

final_layout_t3 <- p_t3_main | layout_b | layout_c
final_layout_t3 <- final_layout_t3 + 
  plot_layout(widths = c(0.8, 3, 2)) + 
  plot_annotation(title = "Task 3: Mu-like Phage Analysis", 
                  subtitle = "B: vFAM (n>10) | C: Hosts from strict vFAMs (n>10 & % in vFAM >= 80%)",
                  theme = theme(plot.title = element_text(size = 18, face = "bold")))

ggsave("Task3_Mu_like_Analysis.pdf", final_layout_t3, width = 18, height = 6)

cat("\nAll customized tasks completed successfully!\n")
