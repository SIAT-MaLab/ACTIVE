
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(foreach)
  library(doParallel)
})

BASE_OUTPUT_DIR <- "host_viral_diversity_analysis_results"
RARE_DEPTH <- 20       
N_ITERATIONS <- 1000   

num_cores <- parallel::detectCores() - 6
if(num_cores < 1) num_cores <- 1
registerDoParallel(cores = num_cores)

dir_mean <- file.path(BASE_OUTPUT_DIR, "1_Mean_Method")
dir_rare <- file.path(BASE_OUTPUT_DIR, "2_Rarefaction_Method")
if (!dir.exists(dir_mean)) dir.create(dir_mean, recursive = TRUE)
if (!dir.exists(dir_rare)) dir.create(dir_rare, recursive = TRUE)
phylum_colors <- c(
  "Actinomycetota"    = "#E41A1C",
  "Bacillota"         = "#377EB8",
  "Bacillota_A"       = "#4DAF4A",
  "Bacillota_B"       = "#984EA3",
  "Bacillota_C"       = "#FF7F00",
  "Bacillota_I"       = "#FFFF33",
  "Bacteroidota"      = "#A65628",
  "Desulfobacterota"  = "#F781BF",
  "Fusobacteriota"    = "#999999",
  "Pseudomonadota"    = "#66C2A5",
  "Verrucomicrobiota" = "#FC8D62"
)

main_file <- "merged_phage_stats_taxonomy.tsv"
if (!file.exists(main_file)) stop("missing merged_phage_stats_taxonomy.tsv")
df <- read.delim(main_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "")
df$Activity_Score <- as.numeric(as.character(df$Activity_Score))

if(file.exists("stats_fisher_vOTU.tsv")){
  sig_votu_df <- read.delim("stats_fisher_vOTU.tsv", header=T, sep="\t", check.names=F)
  list_sig_votus <- sig_votu_df %>% filter(Significance == "Significant") %>% pull(vOTU) %>% unique()
} else { stop("missing stats_fisher_vOTU.tsv") }

if(file.exists("stats_fisher_vFAM.tsv")){
  sig_vfam_df <- read.delim("stats_fisher_vFAM.tsv", header=T, sep="\t", check.names=F)
  list_sig_vfams <- sig_vfam_df %>% filter(Significance == "Significant") %>% pull(vFAM) %>% unique()
} else { stop("missing stats_fisher_vFAM.tsv") }

processed_df <- df %>%
  mutate(
    Is_Active = Activity_Score >= 0.7,
    Is_Sig_Active_vOTU = (vOTU %in% list_sig_votus) & Is_Active,
    Is_Sig_Active_vFAM = (vFAM %in% list_sig_vfams) & Is_Active
  ) %>%
  filter(!is.na(Phylum) & Phylum != "" & !is.na(Family) & Family != "" & !is.na(Genus) & Genus != "")

create_barplot <- function(plot_data, x_col, y_col, fill_col, title, subtitle, y_axis_label) {
  p <- ggplot(plot_data, aes(x = reorder(.data[[x_col]], .data[[y_col]]), y = .data[[y_col]], fill = .data[[fill_col]])) +
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() + 
    scale_fill_manual(values = phylum_colors) +
    labs(title = title, subtitle = subtitle, x = x_col, y = y_axis_label, fill = "Host Phylum") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "italic", size = 10),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )
  return(p)
}

print(">>> [Step 2/3] Mean Diversity...")

run_mean_analysis <- function(data, host_level, metric_col, filter_condition, title_mid, file_prefix) {
  sub_data <- data
  if (filter_condition == "Active") sub_data <- sub_data %>% filter(Is_Active == TRUE)
  if (filter_condition == "Sig_vOTU") sub_data <- sub_data %>% filter(Is_Sig_Active_vOTU == TRUE)
  if (filter_condition == "Sig_vFAM") sub_data <- sub_data %>% filter(Is_Sig_Active_vFAM == TRUE)
  
  stats <- sub_data %>%
    group_by(Genome, across(all_of(c("Phylum", host_level)))) %>%
    summarise(Count_Per_Genome = n_distinct(.data[[metric_col]]), .groups = "drop") %>%
    group_by(across(all_of(c("Phylum", host_level)))) %>%
    summarise(Mean_Count = mean(Count_Per_Genome), N_Genomes = n_distinct(Genome), .groups = "drop") %>%
    filter(N_Genomes >= 5) %>%
    arrange(desc(Mean_Count)) %>%
    slice_head(n = 15)
  
  if(nrow(stats) == 0) return(NULL)
  
  p <- create_barplot(stats, host_level, "Mean_Count", "Phylum", 
                      paste0("Top 15 ", host_level, ": ", title_mid),
                      "Metric: Mean Count per Genome (Min. 5 Genomes)",
                      paste0("Mean Unique ", metric_col))
  ggsave(file.path(dir_mean, paste0(file_prefix, "_", host_level, "_", metric_col, ".pdf")), p, width = 9, height = 7)
}

levels <- c("Family", "Genus", "Species")
for (lvl in levels) {
  run_mean_analysis(processed_df, lvl, "vOTU", "All", "Total vOTU Richness", "A_Total")
  run_mean_analysis(processed_df, lvl, "vGENUS", "All", "Total vGenus Richness", "A_Total")
  run_mean_analysis(processed_df, lvl, "vFAM", "All", "Total vFAM Richness", "A_Total")
  run_mean_analysis(processed_df, lvl, "vOTU", "Active", "Active vOTU Richness", "B_Active")
  run_mean_analysis(processed_df, lvl, "vOTU", "Sig_vOTU", "Significant Active vOTU", "C_Sig")
  run_mean_analysis(processed_df, lvl, "vFAM", "Sig_vFAM", "Significant Active vFAM", "C_Sig")
}


print(">>> [Step 3/3] Parallel Rarefaction...")

run_rarefaction_analysis_parallel <- function(data, host_level, metric_col, filter_condition, title_mid, file_prefix) {
  valid_hosts_df <- data %>%
    group_by(.data[[host_level]], Phylum) %>%
    summarise(n_g = n_distinct(Genome), .groups = "drop") %>%
    filter(n_g >= RARE_DEPTH)
  
  if(nrow(valid_hosts_df) == 0) {
    return(NULL)
  }
  results <- foreach(i = 1:nrow(valid_hosts_df), .combine = rbind, .packages = c("dplyr")) %dopar% {
    
    h_name <- valid_hosts_df[[host_level]][i]
    h_phylum <- valid_hosts_df$Phylum[i]
    sub_df <- data %>% filter(.data[[host_level]] == h_name)
    all_genomes <- unique(sub_df$Genome)
    
    sim_counts <- numeric(N_ITERATIONS)
    for (k in 1:N_ITERATIONS) {
      picked <- sample(all_genomes, RARE_DEPTH, replace = FALSE) 
      
      temp_df <- sub_df %>% filter(Genome %in% picked)
      if (filter_condition == "Active") temp_df <- temp_df %>% filter(Is_Active)
      if (filter_condition == "Sig_vOTU") temp_df <- temp_df %>% filter(Is_Sig_Active_vOTU)
      if (filter_condition == "Sig_vFAM") temp_df <- temp_df %>% filter(Is_Sig_Active_vFAM)
      
      sim_counts[k] <- n_distinct(temp_df[[metric_col]])
    }

    data.frame(
      Host = h_name,
      Phylum = h_phylum,
      Rarefied_Richness = mean(sim_counts),
      stringsAsFactors = FALSE
    )
  }

  if(!is.null(results) && nrow(results) > 0) {
    top_data <- results %>% arrange(desc(Rarefied_Richness)) %>% slice_head(n = 15)
    
    p <- create_barplot(top_data, "Host", "Rarefied_Richness", "Phylum",
                        paste0("Top 15 ", host_level, ": ", title_mid, " (Rarefied)"),
                        paste0("Rarefaction: Depth=", RARE_DEPTH, ", Iterations=", N_ITERATIONS),
                        paste0("Rarefied ", metric_col, " Richness"))
    
    ggsave(file.path(dir_rare, paste0(file_prefix, "_", host_level, "_", metric_col, ".pdf")), p, width = 9, height = 7)
  }
}

for (lvl in levels) {
  run_rarefaction_analysis_parallel(processed_df, lvl, "vOTU", "All", "Total vOTU", "A_Rare_Total")
  run_rarefaction_analysis_parallel(processed_df, lvl, "vGENUS", "All", "Total vGenus", "A_Rare_Total")
  run_rarefaction_analysis_parallel(processed_df, lvl, "vFAM", "All", "Total vFAM", "A_Rare_Total")
  
  run_rarefaction_analysis_parallel(processed_df, lvl, "vOTU", "Active", "Active vOTU", "B_Rare_Active")
  
  run_rarefaction_analysis_parallel(processed_df, lvl, "vOTU", "Sig_vOTU", "Significant Active vOTU", "C_Rare_Sig")
  run_rarefaction_analysis_parallel(processed_df, lvl, "vFAM", "Sig_vFAM", "Significant Active vFAM", "C_Rare_Sig")
}


stopImplicitCluster()

