
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(RColorBrewer)
})

out_dir <- "GLMM_Ultimate_All_In_One_Outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

dir_tables   <- file.path(out_dir, "1_Tables_and_Stats")
dir_variance <- file.path(out_dir, "2_Variance_Partitioning_Plots")
dir_barplots <- file.path(out_dir, "3_Traditional_Barplots_by_Rank")
dir_lollipop <- file.path(out_dir, "4_CrossLevel_Lollipop_Charts")

invisible(lapply(c(dir_tables, dir_variance, dir_barplots, dir_lollipop), function(x) if(!dir.exists(x)) dir.create(x)))
theme_set(theme_bw() + theme(text = element_text(size = 12, color = "black"), axis.text = element_text(color = "black"), panel.grid.minor = element_blank()))
custom_phylum_colors <- c(
  "Bacillota_A"       = "#dce8ef", "Bacillota_C"       = "#add2e5",
  "Bacillota_I"       = "#34666b", "Bacillota"         = "#7abbce",
  "Bacillota_B"       = "#769499", "Actinomycetota"    = "#f5cfa6",
  "Bacteroidota"      = "#efe3ef", "Verrucomicrobiota" = "#d4d68a",
  "Pseudomonadota"    = "#e68f9f", "Desulfobacterota"  = "#d5d5d6",
  "Fusobacteriota"    = "#231815", "Source_Batch"      = "#999999" 
)

df_phage <- read_tsv("../prophage_stats_final.tsv", show_col_types = FALSE)
df_quality <- read_tsv("../allcheckm2_quality.tsv", col_names = FALSE, show_col_types = FALSE)
colnames(df_quality) <- c("Genome_ID", "Completeness", "Contamination", "Model",
                          "Transl_Table", "Coding_Density", "Contig_N50",
                          "Avg_Gene_Len", "Genome_Size", "GC", "Total_CDS", "Total_Contigs", "Max_Contig", "Notes")

df_all <- inner_join(df_phage, df_quality, by = "Genome_ID") %>%
  filter(!is.na(Phylum), !is.na(Genus), !is.na(Species)) %>%
  mutate(
    Prophage_Count = as.numeric(replace_na(Prophage_Count, 0)),
    Active_Prophage_Count = as.numeric(replace_na(Active_Prophage_Count, 0)),
    Is_Lysogen = ifelse(Prophage_Count > 0, 1, 0),
    Log_Prophage_Count = ifelse(Prophage_Count > 0, log(Prophage_Count), NA),
    Is_Active_Carrier = ifelse(Active_Prophage_Count > 0, 1, 0),
    
    across(c(Phylum, Class, Order, Family, Genus, Species, source), as.factor),
    z_Comp = as.numeric(scale(Completeness)),
    z_logN50 = as.numeric(scale(log10(Contig_N50)))
  )

df_lysogens <- df_all %>% filter(Prophage_Count > 0)

create_sample_dict <- function(data_subset) {
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  bind_rows(lapply(ranks, function(r) {
    if (r == "Phylum") {
      res <- data_subset %>% group_by(Phylum) %>% summarise(Sample_Size = n(), .groups = "drop") %>%
        mutate(Taxa_Name = as.character(Phylum), Parent_Phylum = as.character(Phylum))
    } else {
      res <- data_subset %>% group_by(!!sym(r), Phylum) %>% summarise(Sample_Size = n(), .groups = "drop") %>%
        rename(Taxa_Name = !!sym(r), Parent_Phylum = Phylum) %>% mutate(Taxa_Name = as.character(Taxa_Name), Parent_Phylum = as.character(Parent_Phylum))
    }
    res %>% mutate(Rank = r) %>% group_by(Taxa_Name, Rank) %>% slice_max(Sample_Size, n = 1, with_ties = FALSE) %>% ungroup()
  }))
}
sample_dict_prev <- create_sample_dict(df_all)
sample_dict_abun_ind <- create_sample_dict(df_lysogens)

run_glmm_engine <- function(data, formula, family_type, trait_name, sample_dict) {
  if (family_type == "binomial") {
    model <- glmer(formula, data = data, family = binomial, nAGQ = 0, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  } else {
    model <- lmer(formula, data = data, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  }
  vcov_df <- as.data.frame(VarCorr(model)) %>% filter(is.na(var1) | var1 == "(Intercept)") %>%
    mutate(Trait = trait_name, Group = str_split_fixed(grp, ":", 2)[,1], Group = ifelse(Group == "source", "Isolation Source (Batch)", Group),
           Group = ifelse(grp == "Residual", "Residual (Unexplained)", Group), Variance = vcov, Percent_Variance = (Variance / sum(Variance)) * 100) %>%
    select(Trait, Grouping_Factor=grp, Group, Variance, Percent_Variance)
  
  re_df <- as.data.frame(ranef(model)) %>% filter(term == "(Intercept)") %>%
    mutate(Trait = trait_name, Rank = sapply(str_split(grpvar, ":"), `[`, 1), Taxa_Name = sapply(str_split(grp, ":"), `[`, 1), Effect_Size = condval, Std_Error = condsd) %>%
    filter(Rank %in% c("Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% select(Trait, Rank, Taxa_Name, Effect_Size, Std_Error)
  
  final_df <- re_df %>% left_join(sample_dict, by = c("Taxa_Name", "Rank")) %>% mutate(Sample_Size = replace_na(Sample_Size, 1), Parent_Phylum = as.character(replace_na(Parent_Phylum, "Unknown")))
  return(list(variance = vcov_df, effect = final_df))
}

base_cov <- " ~ 1 + z_Comp + z_logN50 + (1|source) + (1|Phylum/Class/Order/Family/Genus/Species)"
res_prev <- run_glmm_engine(df_all, as.formula(paste("Is_Lysogen", base_cov)), "binomial", "Prophage Carriage", sample_dict_prev)
res_abun <- run_glmm_engine(df_lysogens, as.formula(paste("Log_Prophage_Count", base_cov)), "gaussian", "Prophage Burden", sample_dict_abun_ind)
res_ind  <- run_glmm_engine(df_lysogens, as.formula(paste("Is_Active_Carrier", base_cov)), "binomial", "Active Prophage Carriage", sample_dict_abun_ind)


df_var_all <- bind_rows(res_prev$variance, res_abun$variance, res_ind$variance)
write_tsv(df_var_all, file.path(dir_tables, "1_Master_Variance_Partitioning_Stats.tsv"))

rank_order <- c("Residual (Unexplained)", "Isolation Source (Batch)", "Phylum", "Class", "Order", "Family", "Genus", "Species")
df_var_all$Group <- factor(df_var_all$Group, levels = rank_order)
p_var <- ggplot(df_var_all, aes(x = Trait, y = Percent_Variance, fill = Group)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.55) + scale_fill_brewer(palette = "Set3", direction = -1) +
  labs(title = "Variance Partitioning across Key Biological Traits", y = "Percentage of Modeled Variance (%)", x = "") + theme(legend.title = element_blank())
ggsave(file.path(dir_variance, "Variance_Partitioning_Stacked_Plot.pdf"), p_var, width = 8, height = 7)

df_eff_all <- bind_rows(res_prev$effect, res_abun$effect, res_ind$effect)
write_tsv(df_eff_all, file.path(dir_tables, "2_Master_EffectSize_and_BLUPs_Dictionary.tsv"))

plot_classic_bars <- function(df_data, trait_label, target_rank) {
  subset_df <- df_data %>% filter(Trait == trait_label, Rank == target_rank, Sample_Size >= 5) %>% arrange(desc(Effect_Size))
  if(nrow(subset_df) == 0) return(NULL)
  
  if(nrow(subset_df) > 20) {
    plot_df <- bind_rows(head(subset_df, 10) %>% mutate(Status = "Top 10 (Highest Intrinsic Preference)"), tail(subset_df, 10) %>% mutate(Status = "Bottom 10 (Lowest Intrinsic Preference)"))
  } else { plot_df <- subset_df %>% mutate(Status = "All Filtered Taxa") }
  
  plot_df$Taxa_Name <- factor(plot_df$Taxa_Name, levels = plot_df$Taxa_Name[order(plot_df$Effect_Size)])
  plot_df$Status <- factor(plot_df$Status, levels = c("Top 10 (Highest Intrinsic Preference)", "Bottom 10 (Lowest Intrinsic Preference)", "All Filtered Taxa"))
  plot_df <- plot_df %>% mutate(Color_Phylum = ifelse(Parent_Phylum %in% names(custom_phylum_colors), Parent_Phylum, "Source_Batch"))
  safe_trait <- gsub(" \\(.*\\)", "", trait_label) 
  
  p <- ggplot(plot_df, aes(x = Taxa_Name, y = Effect_Size, fill = Color_Phylum)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.75) +
    geom_errorbar(aes(ymin = Effect_Size - Std_Error, ymax = Effect_Size + Std_Error), width=0.2, linewidth=0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth=0.8) + coord_flip() +
    facet_wrap(~Status, scales = "free_y", ncol = 1) + scale_fill_manual(values = custom_phylum_colors) +
    labs(title = paste(safe_trait, "- Absolute Bias-Free Ranking by", target_rank), x = "", y = "Effect Size (Model BLUPs)", fill = "Parent Phylum") +
    theme(strip.text = element_text(face = "bold", size = 11), axis.text.y = element_text(face = "italic"))
  ggsave(file.path(dir_barplots, paste0(safe_trait, "_", target_rank, "_Barplot.pdf")), p, width = 9, height = 8.5)
}
ranks_to_plot <- c("Phylum", "Genus", "Species")
for (t in unique(df_eff_all$Trait)) { for (r in ranks_to_plot) { suppressWarnings(plot_classic_bars(df_eff_all, t, r)) } }



plot_lollipop_cross <- function(plot_df, trait_name, top_n = 10) {
  sdf <- plot_df %>% filter(Trait == trait_name, Sample_Size >= 5) %>% mutate(Label = paste0(Taxa_Name, " (", Rank, ")")) %>% arrange(desc(Effect_Size))
  df_final <- bind_rows(head(sdf, top_n), tail(sdf, top_n)) %>% arrange(Effect_Size) %>% mutate(Label = factor(Label, levels = unique(Label)))
  safe_trait <- gsub(" \\(.*\\)", "", trait_name)
  df_final <- df_final %>% mutate(Color_Phylum = ifelse(Parent_Phylum %in% names(custom_phylum_colors), Parent_Phylum, "Source_Batch"))
  
  p <- ggplot(df_final, aes(x = Effect_Size, y = Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth=0.8) +
    geom_segment(aes(x = 0, xend = Effect_Size, y = Label, yend = Label), color = "grey75", linewidth = 0.6) +
    geom_point(aes(size = Sample_Size, color = Color_Phylum), alpha = 0.9) +
    scale_color_manual(values = custom_phylum_colors) +  # 替换这里为门水平专属调色板
    scale_size_continuous(range = c(2, 12), breaks = c(10, 50, 100, 500, 1000)) +
    labs(title = paste("Top & Bottom Taxa for", safe_trait, "(Cross-Level)"),
         x = "Intrinsic Effect Size (Log Deviation from Grand Mean 0)", y = "",
         size = "Sample Size", color = "Parent Phylum") +
    theme(panel.grid.major.y = element_line(linetype = "dotted", color = "grey85"),
          plot.title = element_text(face = "bold", hjust = 0.5)) + coord_cartesian(clip = "off")
  
  ggsave(file.path(dir_lollipop, paste0(safe_trait, "_CrossLevel_Lollipop.pdf")), p, width = 11.5, height = 9.5)
}

for (t in unique(df_eff_all$Trait)) {
  suppressWarnings(plot_lollipop_cross(df_eff_all, t, top_n = 15))
}
