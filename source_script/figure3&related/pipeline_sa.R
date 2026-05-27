#!/usr/bin/env Rscript

# ==============================================================================
## Function: Advanced Stratified Stats (CMH + Permutation/Exact) Evaluation
#           with Global FDR + Classic Volcano + Pie Chart
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(parallel)     # Native multithreading
})

# ==============================================================================
# 0. Command Line Arguments & Global Configuration
# ==============================================================================
# Parse Command Line Arguments
args <- commandArgs(trailingOnly = TRUE)
PVAL_METHOD <- "permutation"  # Default method
if ("-f" %in% args || "--fisher" %in% args) {
  PVAL_METHOD <- "exact"
}

ACTIVE_THRESHOLD <- 0.7
MIN_GENOME_COUNT <- 5
THREADS <- 20
N_PERMUTATIONS <- 10000

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

cat("========================================================================\n")
cat(sprintf(">>> Initializing Pipeline...\n"))
cat(sprintf(">>> Statistical Method: %s\n", toupper(PVAL_METHOD)))
if(PVAL_METHOD == "permutation") cat(sprintf(">>> Permutations: %d\n", N_PERMUTATIONS))
cat("========================================================================\n")

# ==============================================================================
# 1. Data Loading and Cleaning
# ==============================================================================
cat("\n>>> [1/7] Loading and cleaning data...\n")
raw_meta <- read_tsv("merged_phage_stats_taxonomy.tsv", show_col_types = FALSE)
source_meta <- read_tsv("prophage_stats_final.tsv", show_col_types = FALSE) %>%
  select(Genome_ID, source) %>% distinct()

df_clean <- raw_meta %>%
  left_join(source_meta, by = c("Genome" = "Genome_ID")) %>%
  mutate(Activity_Score = suppressWarnings(as.numeric(as.character(Activity_Score)))) %>%
  mutate(Activity_Score = replace_na(Activity_Score, 0)) %>%
  mutate(Is_Active = ifelse(Activity_Score >= ACTIVE_THRESHOLD, "Yes", "No")) %>%
  mutate(Active_Binary = ifelse(Is_Active == "Yes", 1, 0)) %>%
  mutate(across(c(Species, Genus, Family, Order, Class, Phylum),
                ~ na_if(str_trim(.), "") %>% na_if("s__") %>% na_if("g__") %>% na_if("f__") %>%
                  na_if("Unclassified") %>% na_if("unknown")))

get_mode <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if(length(x) == 0) return(NA)
  names(table(x))[which.max(table(x))]
}

# ==============================================================================
# 2. Advanced Statistical Engine (CMH + Permutation/Exact)
# ==============================================================================
cat(sprintf("\n>>> [2/7] Starting Statistical Engine with %d threads...\n", THREADS))

run_stats_analysis <- function(data, level_name, method = PVAL_METHOD, n_perms = N_PERMUTATIONS) {
  data_lvl <- data %>% filter(!is.na(!!sym(level_name)) & !!sym(level_name) != "")
  
  global_bg <- data_lvl %>%
    group_by(source) %>%
    summarise(Total_Active = sum(Active_Binary), Total_Count = n(), Total_Inactive = Total_Count - Total_Active, .groups = "drop")
  
  taxon_summary <- data_lvl %>%
    group_by(!!sym(level_name)) %>%
    summarise(
      Genome_Count = n(), Count_Active = sum(Active_Binary), Num_Sources = n_distinct(source),
      Max_Intensity = max(Activity_Score), Phylum = get_mode(Phylum), .groups = "drop"
    ) %>% filter(Genome_Count >= MIN_GENOME_COUNT)
  
  taxa_list <- taxon_summary[[level_name]]
  
  results_list <- mclapply(taxa_list, function(taxon) {
    taxon_df <- data_lvl %>% filter(!!sym(level_name) == taxon)
    
    combined <- taxon_df %>%
      group_by(source) %>%
      summarise(a = sum(Active_Binary), c = n() - sum(Active_Binary), .groups = "drop") %>%
      left_join(global_bg, by = "source") %>%
      mutate(b = Total_Active - a, d = Total_Inactive - c) %>%
      filter((a + c) > 0 & (b + d) > 0)
    
    if (nrow(combined) == 0) return(tibble(!!sym(level_name) := taxon, OR = 0.01, CI_Lower = 0.01, CI_Upper = 0.01, p_value = 0.99, Test_Method = "No_Data"))
    
    has_zero <- any(combined$a == 0 | combined$b == 0 | combined$c == 0 | combined$d == 0)
    adj_comb <- combined
    if (has_zero) {
      adj_comb <- adj_comb %>% mutate(a = a + 0.5, b = b + 0.5, c = c + 0.5, d = d + 0.5)
    }
    
    a_adj <- adj_comb$a; b_adj <- adj_comb$b; c_adj <- adj_comb$c; d_adj <- adj_comb$d
    n_adj <- a_adj + b_adj + c_adj + d_adj
    
    or_num <- sum((a_adj * d_adj) / n_adj)
    or_den <- sum((b_adj * c_adj) / n_adj)
    mh_or  <- or_num / or_den

    R <- (a_adj * d_adj) / n_adj; S <- (b_adj * c_adj) / n_adj
    var_logor <- sum((a_adj + d_adj)*R/n_adj^2)/(2*sum(R)^2) + sum((a_adj + d_adj)*S + (c_adj + b_adj)*R)/(2*sum(R)*sum(S)) + sum((c_adj + b_adj)*S/n_adj^2)/(2*sum(S)^2)
    se_logor <- sqrt(var_logor)
    ci_low <- exp(log(mh_or) - 1.96 * se_logor)
    ci_up  <- exp(log(mh_or) + 1.96 * se_logor)
    
    num_strata <- nrow(combined)
    p_val <- NA
    method_used <- method
    
    if (method == "exact") {
      if (num_strata == 1) {
        mat <- matrix(c(combined$a, combined$b, combined$c, combined$d), nrow = 2)
        test_res <- tryCatch(fisher.test(mat), error = function(e) NULL)
        if (!is.null(test_res)) p_val <- test_res$p.value
      } else {
        arr <- array(NA, dim = c(2, 2, num_strata))
        for(i in 1:num_strata) arr[,,i] <- matrix(c(combined$a[i], combined$b[i], combined$c[i], combined$d[i]), nrow = 2)
        test_res <- tryCatch(mantelhaen.test(arr, exact = TRUE), error = function(e) tryCatch(mantelhaen.test(arr, exact = FALSE), error = function(e) NULL))
        if (!is.null(test_res)) p_val <- test_res$p.value
      }
    } else if (method == "permutation") {
      obs_a <- sum(combined$a)
      sim_a_matrix <- matrix(0, nrow = n_perms, ncol = num_strata)
      
      for (i in 1:num_strata) {
        m_i <- combined$a[i] + combined$c[i]
        n_i <- combined$b[i] + combined$d[i]
        k_i <- combined$a[i] + combined$b[i]
        sim_a_matrix[, i] <- rhyper(n_perms, m = m_i, n = n_i, k = k_i)
      }
      
      sim_sum_a <- rowSums(sim_a_matrix)
      mu <- mean(sim_sum_a)
      d_obs <- abs(obs_a - mu)
      
      p_val <- (sum(abs(sim_sum_a - mu) >= d_obs) + 1) / (n_perms + 1)
    }
    
    if (is.na(p_val)) p_val <- 0.99
    
    tibble(!!sym(level_name) := taxon, Test_Method = method_used, OR = mh_or, CI_Lower = ci_low, CI_Upper = ci_up, p_value = p_val)
  }, mc.cores = THREADS)
  
  stats_results <- bind_rows(results_list)
  
  final_df <- taxon_summary %>% left_join(stats_results, by = level_name) %>% filter(!is.na(p_value)) %>%
    mutate(
      Prevalence_Pct = (Count_Active / Genome_Count) * 100, log2_OR = log2(OR),
      p.adj = p.adjust(p_value, method = "fdr"),
      Significance_Type = case_when(p.adj < 0.05 & OR > 1 ~ "Significant", TRUE ~ "Not Significant")
    ) %>% arrange(desc(Significance_Type), p.adj)
  
  return(final_df)
}

stats_votu <- run_stats_analysis(df_clean, "vOTU")
stats_vfam <- run_stats_analysis(df_clean, "vFAM")

write_tsv(stats_votu, sprintf("stats_baseline_vOTU_%s.tsv", toupper(PVAL_METHOD)))
write_tsv(stats_vfam, sprintf("stats_baseline_vFAM_%s.tsv", toupper(PVAL_METHOD)))

cat("\n>>> [3/7] Calculating cross-host capabilities...\n")

votu_cross_stats <- df_clean %>%
  filter(!is.na(vOTU) & vOTU != "") %>%
  group_by(vOTU) %>%
  summarise(
    Member_Count = n(), n_Species = n_distinct(Species, na.rm=TRUE), n_Genus = n_distinct(Genus, na.rm=TRUE),
    n_Family = n_distinct(Family, na.rm=TRUE), n_Order = n_distinct(Order, na.rm=TRUE),
    n_Class = n_distinct(Class, na.rm=TRUE), n_Phylum = n_distinct(Phylum, na.rm=TRUE),
    Max_Activity = max(Activity_Score, na.rm=TRUE), .groups = "drop"
  ) %>%
  mutate(
    Cross_Category = case_when(
      n_Phylum > 1 ~ "Cross-Phylum", n_Class > 1 ~ "Cross-Class", n_Order > 1 ~ "Cross-Order",
      n_Family > 1 ~ "Cross-Family", n_Genus > 1 ~ "Cross-Genus", n_Species > 1 ~ "Cross-Species", TRUE ~ "Non-crossing"
    ),
    Cross_Category = factor(Cross_Category, levels = level_order),
    Is_Generalist = if_else(n_Species > 1, "Generalist", "Specialist")
  )

df_unified <- votu_cross_stats %>%
  left_join(stats_votu %>% select(vOTU, Significance_Type, Phylum), by = "vOTU") %>%
  mutate(
    Significance_Type = replace_na(Significance_Type, "Not Significant"),
    Analysis_Group = case_when(
      Significance_Type == "Significant" ~ "SAvOTU",
      Max_Activity >= ACTIVE_THRESHOLD ~ "Other Active vOTU",
      TRUE ~ "Inactive vOTU"
    ),
    Analysis_Group = factor(Analysis_Group, levels = group_levels)
  )

cat("\n>>> [4/7] Exporting comprehensive detailed tables...\n")
df_full_detailed <- df_unified %>%
  select(vOTU, Phylum, Member_Count, Cross_Category, Is_Generalist, Significance_Type, Analysis_Group, Max_Activity, n_Phylum:n_Species) %>%
  arrange(desc(Cross_Category), desc(Member_Count))
write_csv(df_full_detailed, sprintf("Detailed_Cross_Host_vOTUs_%s.csv", toupper(PVAL_METHOD)))

cat("\n>>> [5/7] Generating Bubble Plots (Prevalence Variant and Classic Volcano)...\n")

prep_bubble_data <- function(stats_data, level_name) {
  if(level_name == "vOTU") {
    plot_df <- stats_data %>%
      left_join(df_full_detailed %>% select(vOTU, Cross_Category), by = "vOTU") %>%
      mutate(Is_Cross_Host = !is.na(Cross_Category) & Cross_Category != "Non-crossing", Cross_Label = case_when(Significance_Type == "Not Significant" ~ "No Cross / Not Sig", Is_Cross_Host ~ "Cross-Host", TRUE ~ "No Cross / Not Sig"))
  } else {
    plot_df <- stats_data %>% mutate(Cross_Label = "No Cross / Not Sig")
  }

  plot_df %>% mutate(
    Color_Group = ifelse(Significance_Type == "Significant", Phylum, "Not Significant"),
    Plot_Order = ifelse(Significance_Type == "Significant", 1, 0),
    log_FDR = -log10(ifelse(p.adj == 0, .Machine$double.xmin, p.adj))
  ) %>% arrange(Plot_Order)
}

draw_prevalence_bubble <- function(stats_data, level_name) {
  plot_df <- prep_bubble_data(stats_data, level_name)
  actual_phyla <- unique(plot_df$Phylum[plot_df$Significance_Type == "Significant"])
  actual_phyla <- actual_phyla[!is.na(actual_phyla)]
  current_palette <- setNames(rep("#bdbdbd", length(actual_phyla)), actual_phyla)
  for(p in actual_phyla) { if(p %in% names(PHYLUM_COLORS)) current_palette[p] <- PHYLUM_COLORS[p] }
  final_colors <- c(current_palette, "Not Significant" = "grey95")
  top_hits <- plot_df %>% filter(Significance_Type == "Significant") %>% arrange(p.adj, desc(log2_OR)) %>% head(10)

  p <- ggplot(plot_df, aes(x = log_FDR, y = Prevalence_Pct)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.6) +
    geom_point(aes(size = Genome_Count, fill = Color_Group, color = Cross_Label), stroke = STROKE_WIDTH, shape = 21, alpha = 0.85) +
    scale_fill_manual(values = final_colors, name = "Host Phylum") +
    scale_color_manual(name = "Cross-Host Status", values = c("Cross-Host" = STROKE_COLOR_CROSS, "No Cross / Not Sig" = STROKE_COLOR_DEFAULT)) +
    scale_size_continuous(range = c(2, 12), name = "Genome Count") +
    labs(title = paste0("Activity Landscape (", toupper(PVAL_METHOD), "): ", level_name), subtitle = paste0("Significantly Active: n=", sum(plot_df$Significance_Type == "Significant")), x = "-log10(FDR)", y = "% of Genome Active (Prevalence)") +
    theme_minimal() + theme(aspect.ratio = 0.75, plot.title = element_text(size=16, face="bold"), panel.border = element_rect(colour="black", fill=NA)) +
    guides(fill = guide_legend(override.aes = list(shape=21, size=5, color=STROKE_COLOR_DEFAULT)))
  if(nrow(top_hits) > 0) p <- p + geom_text_repel(data = top_hits, aes(label = !!sym(level_name)), size=4.5, fontface="bold", min.segment.length=0, segment.color="black", box.padding=0.8, max.overlaps=Inf)
  return(p)
}

draw_classic_volcano <- function(stats_data, level_name) {
  plot_df <- prep_bubble_data(stats_data, level_name)
  actual_phyla <- unique(plot_df$Phylum[plot_df$Significance_Type == "Significant"])
  actual_phyla <- actual_phyla[!is.na(actual_phyla)]
  current_palette <- setNames(rep("#bdbdbd", length(actual_phyla)), actual_phyla)
  for(p in actual_phyla) { if(p %in% names(PHYLUM_COLORS)) current_palette[p] <- PHYLUM_COLORS[p] }
  final_colors <- c(current_palette, "Not Significant" = "grey95")
  top_hits <- plot_df %>% filter(Significance_Type == "Significant") %>% arrange(desc(log2_OR), p.adj) %>% head(10)

  p <- ggplot(plot_df, aes(x = log2_OR, y = log_FDR)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
    geom_point(aes(size = Genome_Count, fill = Color_Group, color = Cross_Label), stroke = STROKE_WIDTH, shape = 21, alpha = 0.85) +
    scale_fill_manual(values = final_colors, name = "Host Phylum") +
    scale_color_manual(name = "Cross-Host Status", values = c("Cross-Host" = STROKE_COLOR_CROSS, "No Cross / Not Sig" = STROKE_COLOR_DEFAULT)) +
    scale_size_continuous(range = c(2, 12), name = "Genome Count") +
    labs(title = paste0("Classic Volcano (", toupper(PVAL_METHOD), "): ", level_name), subtitle = paste0("Significantly Active: n=", sum(plot_df$Significance_Type == "Significant")), x = "Effect Size: log2 (Mantel-Haenszel OR)", y = "-log10(FDR)") +
    theme_minimal() + theme(aspect.ratio = 0.75, plot.title = element_text(size=16, face="bold"), panel.border = element_rect(colour="black", fill=NA)) +
    guides(fill = guide_legend(override.aes = list(shape=21, size=5, color=STROKE_COLOR_DEFAULT)))
  if(nrow(top_hits) > 0) p <- p + geom_text_repel(data = top_hits, aes(label = !!sym(level_name)), size=4.5, fontface="bold", min.segment.length=0, segment.color="black", box.padding=0.8, max.overlaps=Inf)
  return(p)
}

bubble_plot1 <- draw_prevalence_bubble(stats_votu, "vOTU") / draw_prevalence_bubble(stats_vfam, "vFAM") + plot_annotation(tag_levels = 'A')
ggsave(sprintf("Combined_Activity_Landscape_%s.pdf", toupper(PVAL_METHOD)), bubble_plot1, width = 12, height = 18, bg = "white")

bubble_plot2 <- draw_classic_volcano(stats_votu, "vOTU") / draw_classic_volcano(stats_vfam, "vFAM") + plot_annotation(tag_levels = 'A')
ggsave(sprintf("Combined_Classic_Volcano_%s.pdf", toupper(PVAL_METHOD)), bubble_plot2, width = 12, height = 18, bg = "white")
# --- Plot Type C: X = log2_OR, Y = Prevalence_Pct ---
draw_prevalence_vs_or_bubble <- function(stats_data, level_name) {
  plot_df <- prep_bubble_data(stats_data, level_name)
  actual_phyla <- unique(plot_df$Phylum[plot_df$Significance_Type == "Significant"])
  actual_phyla <- actual_phyla[!is.na(actual_phyla)]
  current_palette <- setNames(rep("#bdbdbd", length(actual_phyla)), actual_phyla)
  for(p in actual_phyla) { if(p %in% names(PHYLUM_COLORS)) current_palette[p] <- PHYLUM_COLORS[p] }
  final_colors <- c(current_palette, "Not Significant" = "grey95")

  top_hits <- plot_df %>% filter(Significance_Type == "Significant") %>% arrange(desc(log2_OR), desc(Prevalence_Pct)) %>% head(10)

  p <- ggplot(plot_df, aes(x = log2_OR, y = Prevalence_Pct)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.6) +
    geom_point(aes(size = Genome_Count, fill = Color_Group, color = Cross_Label), stroke = STROKE_WIDTH, shape = 21, alpha = 0.85) +
    scale_fill_manual(values = final_colors, name = "Host Phylum") +
    scale_color_manual(name = "Cross-Host Status", values = c("Cross-Host" = STROKE_COLOR_CROSS, "No Cross / Not Sig" = STROKE_COLOR_DEFAULT)) +
    scale_size_continuous(range = c(2, 12), name = "Genome Count") +
    labs(title = paste0("Prevalence vs Effect Size (", toupper(PVAL_METHOD), "): ", level_name), 
         subtitle = paste0("Significantly Active: n=", sum(plot_df$Significance_Type == "Significant")), 
         x = "Effect Size: log2 (Mantel-Haenszel OR)", 
         y = "% of Genome Active (Prevalence)") +
    theme_minimal() + theme(aspect.ratio = 0.75, plot.title = element_text(size=16, face="bold"), panel.border = element_rect(colour="black", fill=NA)) +
    guides(fill = guide_legend(override.aes = list(shape=21, size=5, color=STROKE_COLOR_DEFAULT)))
    
  if(nrow(top_hits) > 0) p <- p + geom_text_repel(data = top_hits, aes(label = !!sym(level_name)), size=4.5, fontface="bold", min.segment.length=0, segment.color="black", box.padding=0.8, max.overlaps=Inf)
  return(p)
}

bubble_plot3 <- draw_prevalence_vs_or_bubble(stats_votu, "vOTU") / draw_prevalence_vs_or_bubble(stats_vfam, "vFAM") + plot_annotation(tag_levels = 'A')
ggsave(sprintf("Combined_Prevalence_vs_OR_%s.pdf", toupper(PVAL_METHOD)), bubble_plot3, width = 12, height = 18, bg = "white")

cat("\n>>> [6/7] Generating Comprehensive 2-in-1 Analysis...\n")

df_A <- df_unified %>% filter(Analysis_Group == "SAvOTU") %>% count(Cross_Category) %>% mutate(Prop = n/sum(n))
if(nrow(df_A) > 0) {
  pA <- ggplot(df_A, aes(x = Cross_Category, y = n, fill = Cross_Category)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") + geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, Prop*100)), vjust = -0.2, size = 3.5) +
    scale_fill_manual(values = cross_colors) + scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "A. Cross-Host Distribution of SAvOTUs", x = "", y = "Number of SAvOTUs") + theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1, face="bold"))
} else {
  pA <- ggplot() + theme_void() + labs(title = "A. No SAvOTUs Found")
}

df_B <- df_unified %>% group_by(Analysis_Group) %>% summarise(Total_vOTU = n(), Generalists = sum(Is_Generalist == "Generalist"), Prop_Generalist = Generalists/Total_vOTU, .groups="drop")

calc_fisher <- function(g1, g2) {
  dat <- df_unified %>% filter(Analysis_Group %in% c(g1, g2))
  if(length(unique(dat$Analysis_Group)) < 2) return(1)
  fisher.test(table(dat$Analysis_Group, dat$Is_Generalist))$p.value
}

raw_pval_S_O <- calc_fisher("SAvOTU", "Other Active vOTU")
raw_pval_S_I <- calc_fisher("SAvOTU", "Inactive vOTU")
adj_pvals <- p.adjust(c(raw_pval_S_O, raw_pval_S_I), method = "fdr")
pval_S_O <- adj_pvals[1]
pval_S_I <- adj_pvals[2]

format_pval <- function(p) { if(p < 0.0001) "****" else if(p < 0.001) "***" else if(p < 0.01) "**" else if(p < 0.05) "*" else "ns" }
ymax <- max(df_B$Prop_Generalist, na.rm = TRUE)

if(ymax > 0 && "SAvOTU" %in% df_B$Analysis_Group) {
  pB <- ggplot(df_B, aes(x = Analysis_Group, y = Prop_Generalist, fill = Analysis_Group)) +
    geom_bar(stat="identity", width=0.6, color="black") + geom_text(aes(label = percent(Prop_Generalist, accuracy=0.1)), vjust=-0.5, size=4) +
    geom_segment(aes(x=2, xend=3, y=ymax*1.15, yend=ymax*1.15)) + geom_text(aes(x=2.5, y=ymax*1.18), label=sprintf("%s\n(p=%.2e)", format_pval(pval_S_O), pval_S_O), size=3) +
    geom_segment(aes(x=1, xend=3, y=ymax*1.35, yend=ymax*1.35)) + geom_text(aes(x=2, y=ymax*1.38), label=sprintf("%s\n(p=%.2e)", format_pval(pval_S_I), pval_S_I), size=3) +
    scale_fill_manual(values = group_colors) + scale_y_continuous(labels = percent, limits = c(0, ymax*1.5)) + labs(title = "B. Proportion of Generalists", x="", y="Generalist Proportion (%)") + theme_classic() + theme(legend.position="none", axis.text.x = element_text(angle=15, hjust=1, face="bold"))
} else {
  pB <- ggplot(df_B, aes(x = Analysis_Group, y = Prop_Generalist, fill = Analysis_Group)) + geom_bar(stat="identity", width=0.6, color="black") + theme_classic() + labs(title = "B. Proportion of Generalists (Insufficient Data)")
}

comprehensive_plot <- pA / pB + plot_annotation(title = sprintf("SAvOTU Analysis (%s)", toupper(PVAL_METHOD)))
ggsave(sprintf("Fig_Comprehensive_2Groups_%s.pdf", toupper(PVAL_METHOD)), comprehensive_plot, width = 8, height = 9)

cat("\n>>> [7/7] Generating Phylum Composition Pie Chart for SAvOTUs...\n")

df_pie <- df_unified %>%
  filter(Analysis_Group == "SAvOTU") %>% filter(!is.na(Phylum)) %>%
  group_by(Phylum) %>% summarise(Count = n(), .groups = "drop") %>% arrange(desc(Count)) %>%
  mutate(Prop = Count / sum(Count), Label_Text = ifelse(Prop >= 0.02, paste0(Phylum, "\n", Count, " (", percent(Prop, accuracy = 0.1), ")"), ""))

if(nrow(df_pie) > 0) {
  pie_plot <- ggplot(df_pie, aes(x = "", y = Count, fill = Phylum)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = PHYLUM_COLORS) +
    theme_void() +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10))) +
    labs(title = sprintf("Phylum Composition of SAvOTUs (%s)", toupper(PVAL_METHOD)), fill = "Host Phylum") +
    geom_text(aes(label = Label_Text), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold")
  
  ggsave(sprintf("PieChart_SAvOTU_Phylum_Composition_%s.pdf", toupper(PVAL_METHOD)), pie_plot, width = 8, height = 6, bg = "white")
} else {
  cat("\n[!] No SAvOTUs found. Skipping pie chart generation.\n")
}

cat("\n>>> Pipeline successfully completed! All tables and plots generated.\n")
