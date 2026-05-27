

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lme4)
  library(doParallel)
  library(foreach)
  library(patchwork)
  library(scales)
  library(stringr)
})


ANNOT_FILE    <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_annot/allpharokka_annot.tsv"
METADATA_FILE <- "/home/huangyan/experiment/SPA/3/gut_isolate/prophage/analysis/3_prophage_stats/allpp_checkv_quality/hqphage_metadata.tsv"

OUT_TABLE     <- "Pharokka_Annot_TwoStep_GLM_GLMM_Stats.tsv"
OUT_PLOT      <- "Pharokka_Enrichment_GLMM_TopBottom_Activity.pdf"

GLM_FDR_CUTOFF   <- 0.05  
GLMM_FDR_CUTOFF  <- 0.05  
MIN_PROP          <- 0.01  
N_THREADS         <- 25

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
    vFAM = ifelse(is.na(vFAM) | vFAM == "", "Unknown_vFAM", vFAM),
    location_type = as.factor(Phage_Location_Type)
  )

top_phyla <- get_top_taxa(df_meta$Phylum, n=15)
df_meta <- df_meta %>%
  mutate(
    host_taxa_reg = ifelse(Phylum %in% top_phyla, Phylum, "Other"),
    host_taxa_reg = as.factor(host_taxa_reg)
  ) %>%
  filter(!is.na(log_len) & !is.na(host_taxa_reg))


n_total_phages  <- nrow(df_meta)
n_induced_total <- sum(df_meta$Is_Induced == 1)

cat(sprintf("   - Total modeling phages: %d, Induced: %d\n", n_total_phages, n_induced_total))


df_annot_raw <- fread(ANNOT_FILE, select = c("contig", "annot")) %>%
  rename(Phage_Name = contig) %>%
  mutate(annot_clean = tolower(trimws(annot))) %>%
  filter(!grepl("hypothetical|unknown function", annot_clean)) %>%
  filter(annot_clean != "") %>%
  distinct(Phage_Name, annot_clean) 

annot_counts <- df_annot_raw %>%
  inner_join(df_meta %>% select(Phage_Name), by = "Phage_Name") %>% 
  distinct(Phage_Name, annot_clean) %>%
  count(annot_clean, name = "Total")

valid_annots <- annot_counts %>%
  filter(Total > n_total_phages * MIN_PROP) %>% 
  pull(annot_clean)

cat(sprintf("   - Total features after >1%% proportion filtering: %d\n", length(valid_annots)))

cat(sprintf("Step 2: Parallel GLM Screening using %d threads...\n", N_THREADS))

fixed_effects <- "Has_Annot + host_taxa_reg + location_type + log_len"

cl <- makeCluster(N_THREADS)
registerDoParallel(cl)

glm_results <- foreach(
  i = seq_along(valid_annots),
  .combine = rbind,
  .packages = c("stats", "dplyr", "stringr"),
  .export = c("valid_annots", "df_annot_raw", "df_meta", "fixed_effects"),
  .errorhandling = "remove"
) %dopar% {

  target_annot <- valid_annots[i]
  phages_with_annot <- df_annot_raw %>% filter(annot_clean == target_annot) %>% pull(Phage_Name)

  df_model <- df_meta %>%
    mutate(Has_Annot = ifelse(Phage_Name %in% phages_with_annot, 1, 0))

  total_annot <- sum(df_model$Has_Annot == 1)
  induced_annot <- sum(df_model$Has_Annot == 1 & df_model$Is_Induced == 1)

  res <- tryCatch({
    form <- as.formula(paste0("Is_Induced ~ ", fixed_effects))
    m <- glm(form, data = df_model, family = binomial)

    coefs <- summary(m)$coefficients
    if ("Has_Annot" %in% rownames(coefs)) {
      c(coefs["Has_Annot", "Estimate"], coefs["Has_Annot", "Pr(>|z|)"])
    } else { c(NA, NA) }
  }, error = function(e) { c(NA, NA) })

  data.frame(
    annot_clean = target_annot,
    Total = total_annot,
    Induced = induced_annot,
    GLM_Est_raw = res[1],
    GLM_P = res[2]
  )
}
stopCluster(cl)
glm_results <- glm_results %>%
  mutate(
    GLM_Log2Est = GLM_Est_raw / log(2),
    GLM_FDR = p.adjust(GLM_P, method = "BH"),
    Is_Candidate = !is.na(GLM_FDR) & GLM_FDR < GLM_FDR_CUTOFF
  ) %>%
  select(-GLM_Est_raw)

candidates <- glm_results %>% filter(Is_Candidate) %>% pull(annot_clean)
cat(sprintf("   - GLM Candidates for GLMM validation: %d\n", length(candidates)))

if(length(candidates) > 0) {
  cat(sprintf("Step 3: Parallel GLMM Validation using %d threads...\n", N_THREADS))

  cl <- makeCluster(N_THREADS)
  registerDoParallel(cl)

  glmm_results <- foreach(
    i = seq_along(candidates),
    .combine = rbind,
    .packages = c("lme4", "dplyr", "stringr"),
    .export = c("candidates", "df_annot_raw", "df_meta", "fixed_effects"),
    .errorhandling = "remove"
  ) %dopar% {

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
        c(coefs["Has_Annot", "Estimate"], coefs["Has_Annot", "Pr(>|z|)"])
      } else { c(NA, NA) }
    }, error = function(e) { c(NA, NA) })

    data.frame(
      annot_clean = target_annot,
      GLMM_Est_raw = res[1],
      GLMM_P = res[2]
    )
  }

  stopCluster(cl)

  glmm_results <- glmm_results %>%
    mutate(
      GLMM_Log2Est = GLMM_Est_raw / log(2),
      GLMM_FDR = p.adjust(GLMM_P, method = "BH"),
      Is_Validated = !is.na(GLMM_FDR) & GLMM_FDR < GLMM_FDR_CUTOFF
    ) %>%
    select(-GLMM_Est_raw)

  final_data <- glm_results %>%
    inner_join(glmm_results, by = "annot_clean") %>%
    mutate(
      Cryptic = Total - Induced,
      Rate_Encoded = Induced / Total,
      Rate_Bg = (n_induced_total - Induced) / (n_total_phages - Total),
      Group = case_when(
        Is_Validated & GLMM_Log2Est > 0 ~ "Induced",
        Is_Validated & GLMM_Log2Est < 0 ~ "Cryptic",
        TRUE ~ "NS"
      )
    )

  write_tsv(final_data, OUT_TABLE)
  cat(sprintf("   - Results saved to %s\n", OUT_TABLE))

  cat("Step 4: Generating TOP/Bottom Plot...\n")

  validated_data <- final_data %>% filter(Is_Validated)

  if(nrow(validated_data) > 0) {

    top10_induced <- validated_data %>%
      filter(GLMM_Log2Est > 0) %>%
      arrange(desc(GLMM_Log2Est)) %>%
      slice_head(n = 10)

    bottom10_cryptic <- validated_data %>%
      filter(GLMM_Log2Est < 0) %>%
      arrange(GLMM_Log2Est) %>%
      slice_head(n = 10)

    plot_data <- bind_rows(top10_induced, bottom10_cryptic) %>%
      mutate(
        label_display = str_to_sentence(annot_clean),
        label_display = str_trunc(label_display, 40),
        label_display = fct_reorder(label_display, GLMM_Log2Est),
        Direction = ifelse(GLMM_Log2Est > 0, "Induced-enriched", "Cryptic-enriched")
      )

    rate_long <- plot_data %>%
      select(label_display, Rate_Encoded, Rate_Bg) %>%
      pivot_longer(cols = c(Rate_Encoded, Rate_Bg), names_to = "Type", values_to = "Rate") %>%
      mutate(Type = ifelse(Type == "Rate_Encoded", "Encoded", "Background"))

    p_est <- ggplot(plot_data, aes(x = GLMM_Log2Est, y = label_display)) +
      geom_col(aes(fill = Direction), width = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
      scale_fill_manual(values = c("Induced-enriched" = "#e74c3c", "Cryptic-enriched" = "#3498db")) +
      labs(x = "GLMM Log2 Estimate\n(Adj. for Host/Loc/Len/vFAM)", y = NULL, title = "Effect Size (GLMM)") +
      theme_bw() +
      theme(panel.grid.major.y = element_blank(), legend.position = "none", axis.text.y = element_text(size = 10))

    ggsave(OUT_PLOT, p_est, width = 8, height = 6)
    cat(sprintf("   - Plot saved to %s\n", OUT_PLOT))
  } else {
    cat("   - No validated candidates found. Skipping plot.\n")
  }
} else {
  cat("   - No candidates passed GLM screening. Skipping GLMM and Plot.\n")
}

cat("Done!\n")
