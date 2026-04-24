
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lme4)
  library(patchwork)
  library(scales)
  library(stringr) 
})

ANNOT_FILE    <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_annot/allpharokka_annot.tsv"
METADATA_FILE <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_checkv_quality/hqphage_metadata.tsv"

OUT_TABLE_FULL <- "Pharokka_Annot_GLMM_Stats_Covariates.tsv"
OUT_PLOT_ALL   <- "Pharokka_Enrichment_Plot_A_Comparison_Top10.pdf"
OUT_PLOT_STRICT<- "Pharokka_Enrichment_Plot_B_Strict_Validated_Top25_Covariates.pdf" 

FISHER_FDR_CUTOFF <- 0.05
FISHER_LOG2OR_CUT <- 0
GLMM_P_CUTOFF     <- 0.05
MIN_COUNT         <- 15

cat("Step 1: Loading and Preprocessing Data...\n")


df_meta_raw <- fread(METADATA_FILE)


get_top_taxa <- function(taxa_vec, n=15) {
  counts <- sort(table(taxa_vec), decreasing = TRUE)
  top_names <- names(counts)[1:min(length(counts), n)]
  return(top_names)
}

df_meta <- df_meta_raw %>%
  filter(grepl("Caudoviricetes", ictv_Class)) %>%

  mutate(
    coords = str_extract(Phage_Name, "[0-9]+-[0-9]+$"),
    start_pos = as.numeric(str_extract(coords, "^[0-9]+")),
    end_pos   = as.numeric(str_extract(coords, "[0-9]+$")),
    genome_len = end_pos - start_pos + 1,
    log_len = log10(genome_len) 
  ) %>%

  mutate(
    Is_Induced = ifelse(Activity_Score >= 0.7 & !is.na(Activity_Score), 1, 0),
    Status_Label = ifelse(Is_Induced == 1, "Induced", "Cryptic"),
    vFAM = ifelse(is.na(vFAM) | vFAM == "", "Unknown_vFAM", vFAM),
    location_type = as.factor(Phage_Location_Type) 
  )


top_phyla <- get_top_taxa(df_meta$Phylum, n=15)
df_meta <- df_meta %>%
  mutate(
    host_taxa_reg = ifelse(Phylum %in% top_phyla, Phylum, "Other"),
    host_taxa_reg = as.factor(host_taxa_reg)
  )


df_meta <- df_meta %>% filter(!is.na(log_len) & !is.na(host_taxa_reg))


n_total_phages  <- nrow(df_meta)
n_induced_total <- sum(df_meta$Is_Induced == 1)
n_cryptic_total <- sum(df_meta$Is_Induced == 0)

cat(sprintf("   - Total: %d, Induced: %d, Cryptic: %d\n", n_total_phages, n_induced_total, n_cryptic_total))
cat("   - Covariates added: host_taxa_reg, location_type, log10(genome_len)\n")


df_annot_raw <- fread(ANNOT_FILE, select = c("contig", "annot")) %>%
  rename(Phage_Name = contig) %>%
  mutate(annot_clean = tolower(trimws(annot))) %>%
  filter(!grepl("hypothetical|unknown function", annot_clean)) %>%
  filter(annot_clean != "") %>%
  distinct(Phage_Name, annot_clean)


df_counts <- df_annot_raw %>%
  inner_join(df_meta, by = "Phage_Name") %>%
  count(annot_clean, Status_Label) %>%
  pivot_wider(names_from = Status_Label, values_from = n, values_fill = 0)

if(!"Induced" %in% names(df_counts)) df_counts$Induced <- 0
if(!"Cryptic" %in% names(df_counts)) df_counts$Cryptic <- 0

df_filtered <- df_counts %>%
  mutate(Total = Induced + Cryptic) %>%
  filter(Total >= MIN_COUNT)

cat("Step 2: Fisher Screening...\n")

run_fisher <- function(idx) {
  row <- df_filtered[idx, ]
  mat <- matrix(c(row$Induced, row$Cryptic,
                  n_induced_total - row$Induced, n_cryptic_total - row$Cryptic), nrow=2)
  ft <- fisher.test(mat, conf.level = 0.95)
  
  est_val <- ft$estimate
  if(est_val == 0) est_val <- 1e-10
  if(is.infinite(est_val)) est_val <- 1e10
  
  log2or <- log2(est_val)
  ci_lower <- log2(ft$conf.int[1])
  ci_upper <- log2(ft$conf.int[2])
  
  if(is.infinite(ci_lower)) ci_lower <- -10
  if(is.infinite(ci_upper)) ci_upper <- 10
  
  return(data.frame(
    annot_clean = row$annot_clean,
    Fisher_P = ft$p.value,
    Log2OR = log2or,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper
  ))
}

fisher_res <- do.call(rbind, lapply(1:nrow(df_filtered), run_fisher))

df_fisher_final <- df_filtered %>%
  left_join(fisher_res, by = "annot_clean") %>%
  mutate(
    Fisher_FDR = p.adjust(Fisher_P, method = "BH"),
    Is_Candidate = Fisher_FDR < FISHER_FDR_CUTOFF & abs(Log2OR) > FISHER_LOG2OR_CUT,
    Group = case_when(
      Is_Candidate & Log2OR > 0 ~ "Induced",
      Is_Candidate & Log2OR < 0 ~ "Cryptic",
      TRUE ~ "NS"
    )
  )

candidates <- df_fisher_final %>% filter(Is_Candidate) %>% pull(annot_clean)
cat(sprintf("   - Candidates for GLMM: %d\n", length(candidates)))

cat("Step 3: GLMM Validation (Adjusted for Host, Location, Length)...\n")

glmm_results <- list()


#  Is_Induced ~ Has_Annot + host_taxa_reg + location_type + log_len + (1 | vFAM)
fixed_effects <- "Has_Annot + host_taxa_reg + location_type + log_len"

for (i in seq_along(candidates)) {
  if(i %% 10 == 0) cat(sprintf("Processing %d / %d ...\n", i, length(candidates)))
  
  target_annot <- candidates[i]
  
  phages_with_annot <- df_annot_raw %>% filter(annot_clean == target_annot) %>% pull(Phage_Name)
  
  df_model <- df_meta %>% 
    mutate(Has_Annot = ifelse(Phage_Name %in% phages_with_annot, 1, 0))
  
  res <- tryCatch({
    form <- as.formula(paste0("Is_Induced ~ ", fixed_effects, " + (1 | vFAM)"))

    model <- glmer(form,
                   data = df_model, family = binomial,
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=1e5)), 
                   nAGQ=0)
    
    coefs <- summary(model)$coefficients

    if ("Has_Annot" %in% rownames(coefs)) {
      est <- coefs["Has_Annot", "Estimate"]
      p_val <- coefs["Has_Annot", "Pr(>|z|)"]
      c(est, p_val)
    } else { c(NA, NA) }
  }, error = function(e) { 
    # message(paste("Error in", target_annot, ":", e$message)) 
    c(NA, NA) 
  })
  
  glmm_results[[i]] <- data.frame(
    annot_clean = target_annot,
    GLMM_Log2Est = res[1] / log(2), 
    GLMM_P = res[2]
  )
}

df_glmm <- do.call(rbind, glmm_results)

final_data <- df_fisher_final %>%
  left_join(df_glmm, by = "annot_clean") %>%
  mutate(
    Is_Validated = !is.na(GLMM_P) & GLMM_P < GLMM_P_CUTOFF,
    NegLogP = -log10(GLMM_P),
    NegLogP = ifelse(is.infinite(NegLogP) | NegLogP > 10, 10, NegLogP),
    NegLogP = ifelse(is.na(NegLogP), 0, NegLogP),
    
    # Activity Rate
    Rate_Encoded = Induced / Total,
    Rate_Bg = (n_induced_total - Induced) / (n_total_phages - Total)
  )

write_tsv(final_data, OUT_TABLE_FULL)
cat(sprintf("   - Results saved to %s\n", OUT_TABLE_FULL))


create_dual_plot <- function(data_subset, title_text, subtitle_text) {
  plot_df <- data_subset %>%
    mutate(
      label_display = str_to_sentence(annot_clean),
      label_display = str_trunc(label_display, 40),
      label_display = fct_reorder(label_display, Log2OR)
    )
  
  p_left <- ggplot(plot_df, aes(x = Log2OR, y = label_display)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper, color = Group), height = 0.2) +
    geom_point(aes(color = Group), size = 3) +
    scale_color_manual(values = c("Induced" = "#e74c3c", "Cryptic" = "#3498db")) +
    labs(title = "Prevalence (Fisher)", x = "Log2 Odds Ratio (95% CI)", y = NULL) +
    theme_bw() + theme(panel.grid.major.y = element_blank(), legend.position = "none")
  
  p_right <- ggplot(plot_df, aes(x = NegLogP, y = label_display)) +
    geom_col(aes(fill = Is_Validated), width = 0.6) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("TRUE" = "#2ecc71", "FALSE" = "gray70")) +
    labs(title = "Validation (GLMM)", x = "-log10(P-value)", y = NULL) +
    theme_bw() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
  
  return(p_left + p_right + plot_layout(widths = c(1.2, 0.8)) + plot_annotation(title = title_text, subtitle = subtitle_text))
}


create_tri_plot <- function(data_subset, title_text, subtitle_text) {
  plot_df <- data_subset %>%
    mutate(
      label_display = str_to_sentence(annot_clean),
      label_display = str_trunc(label_display, 40),
      label_display = fct_reorder(label_display, Log2OR)
    )
  
  rate_long <- plot_df %>%
    select(label_display, Rate_Encoded, Rate_Bg) %>%
    pivot_longer(cols = c(Rate_Encoded, Rate_Bg), names_to = "Type", values_to = "Rate") %>%
    mutate(Type = ifelse(Type == "Rate_Encoded", "Encoded", "Background"))
  
  p_left <- ggplot(plot_df, aes(x = Log2OR, y = label_display)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper, color = Group), height = 0.2) +
    geom_point(aes(color = Group), size = 3) +
    scale_color_manual(values = c("Induced" = "#e74c3c", "Cryptic" = "#3498db")) +
    labs(title = "Prevalence (Fisher)", x = "Log2 Odds Ratio (95% CI)", y = NULL) +
    theme_bw() + theme(panel.grid.major.y = element_blank(), legend.position = "none")
  
  p_mid <- ggplot(plot_df, aes(x = NegLogP, y = label_display)) +
    geom_col(aes(fill = Is_Validated), width = 0.6) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("TRUE" = "#2ecc71", "FALSE" = "gray70")) +
    labs(title = "Validation (GLMM)", x = "-log10(P-value)", y = NULL) +
    theme_bw() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
  
  p_right <- ggplot(rate_long, aes(x = Rate, y = label_display)) +
    geom_line(aes(group = label_display), color = "gray80", size = 0.8) +
    geom_point(aes(color = Type), size = 3) +
    scale_color_manual(values = c("Encoded" = "#e74c3c", "Background" = "gray40"), name = "Phage Activity") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Activity Rate", x = "% Induced", y = NULL) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
  
  return((p_left + p_mid + p_right) + plot_layout(widths = c(1.2, 0.6, 1.0)) + plot_annotation(title = title_text, subtitle = subtitle_text))
}


cat("Step 4: Generating Plot A (Top 10 Comparison)...\n")
data_plot_A <- final_data %>%
  filter(Is_Candidate) %>%
  group_by(Group) %>%
  arrange(desc(abs(Log2OR))) %>%
  slice_head(n = 10) %>%
  ungroup()

plot_A <- create_dual_plot(data_plot_A, "Plot A: Comparison (Top 10)", "Fisher candidates regardless of validation")
ggsave(OUT_PLOT_ALL, plot_A, width = 12, height = 8)

cat("Step 5: Generating Plot B (Top 25 Strict + Activity)...\n")
data_plot_B <- final_data %>%
  filter(Is_Candidate & Is_Validated) %>%
  group_by(Group) %>%
  arrange(desc(abs(Log2OR))) %>%
  slice_head(n = 25) %>%
  ungroup()

plot_B <- create_tri_plot(data_plot_B, "Plot B: Validated Features (Top 25)", "Left: Prevalence, Mid: GLMM P-value (Adj. for Host/Loc/Len), Right: Activity")
ggsave(OUT_PLOT_STRICT, plot_B, width = 16, height = 12)

cat("Done!\n")
