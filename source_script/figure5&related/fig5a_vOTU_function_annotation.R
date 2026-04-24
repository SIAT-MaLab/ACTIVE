
library(tidyverse)
library(data.table)
library(lme4)     
library(lmerTest) 

META_FILE    <- "hqphage_metadata.tsv"
PC_MAP_FILE  <- "allpp_PC_map.tsv"
ANNOT_FILE   <- "mmseqs_repseq_annot/final/annotations-Protein_mmseqs_out_rep_seq/final_annotation_summary.tsv"
OUT_PREFIX   <- "Final_TwoStep_"
meta    <- fread(META_FILE)
pc_map  <- fread(PC_MAP_FILE)
annot   <- fread(ANNOT_FILE)

get_mode <- function(x) {
  ux <- unique(na.omit(x))
  if(length(ux) == 0) return("Unknown")
  ux[which.max(tabulate(match(x, ux)))]
}

meta_clean <- meta %>%
  filter(grepl("Caudoviricetes", ictv_Class)) %>%
  mutate(is_active = if_else(!is.na(Activity_Score) & Activity_Score >= 0.7, 1, 0))

otu_stats <- meta_clean %>%
  group_by(vOTU) %>%
  summarise(
    total_genomes = n(),
    active_count = sum(is_active),
    inactive_count = n() - sum(is_active),
    host_phylum = get_mode(Phylum),
    host_class  = get_mode(Class),
    vFAM        = get_mode(vFAM),
    vGENUS      = get_mode(vGENUS),
    location_type = get_mode(Phage_Location_Type)
  ) %>%
  ungroup() %>%
  filter(total_genomes >= 5) 

target_taxa_col <- "host_phylum"
top_taxa <- names(sort(table(otu_stats[[target_taxa_col]]), decreasing = TRUE))[1:15]

otu_stats <- otu_stats %>%
  mutate(host_taxa_reg = if_else(.data[[target_taxa_col]] %in% top_taxa, .data[[target_taxa_col]], "Other")) %>%
  mutate(across(c(host_taxa_reg, location_type, vFAM, vGENUS), as.factor))

RANDOM_EFF <- "vFAM"

besthit_dict <- annot %>%
  filter(best_hit != "") %>%
  group_by(best_hit) %>% slice(1) %>% ungroup() %>%
  select(Feature = best_hit, Product = product)

pc_map_linked <- pc_map %>%
  mutate(Phage_Name = str_remove(protein_id, "_CDS_[0-9]+$")) %>%
  inner_join(meta_clean %>% select(Phage_Name, vOTU), by = "Phage_Name")

get_prev_matrix <- function(link_df, otu_df) {
  link_df %>%
    distinct(vOTU, Feature, Phage_Name) %>%
    count(vOTU, Feature, name = "hits") %>%
    inner_join(otu_df %>% select(vOTU, total_genomes), by = "vOTU") %>%
    mutate(prev = hits / total_genomes) %>%
    mutate(prev = if_else(prev > 1, 1, prev)) %>%

    group_by(Feature) %>%
    filter(sum(prev >= 0.4) >= 5) %>%
    ungroup() %>%
    select(vOTU, Feature, prev) %>%
    pivot_wider(names_from = Feature, values_from = prev, values_fill = 0)
}


run_two_step_analysis <- function(feature_mat, otu_df, feat_type) {
  df <- otu_df %>% inner_join(feature_mat, by = "vOTU")
  features <- setdiff(names(feature_mat), "vOTU")
  fixed_f <- "host_taxa_reg + location_type + log10(total_genomes)"
  glm_res <- data.frame()
  for(i in seq_along(features)) {
    f <- features[i]
    if(i %% 200 == 0) message(paste0("  GLM Progress: ", i, "/", length(features)))
    
    form <- paste0("cbind(active_count, inactive_count) ~ `", f, "` + ", fixed_f)
    try({
      m <- glm(as.formula(form), data = df, family = quasibinomial)
      coef <- summary(m)$coefficients
      glm_res <- rbind(glm_res, data.frame(
        Feature = f, 
        GLM_Est = coef[2, 1], 
        GLM_P = coef[2, 4]
      ))
    }, silent = TRUE)
  }
  if(nrow(glm_res) == 0) return(NULL)
  glm_res$GLM_FDR <- p.adjust(glm_res$GLM_P, method = "BH")
  candidates <- glm_res %>% filter(GLM_FDR <= 0.05) %>% pull(Feature)
  if(length(candidates) == 0) return(NULL)
  glmm_res <- data.frame()
  for(i in seq_along(candidates)) {
    f <- candidates[i]
    if(i %% 10 == 0) message(paste0("  GLMM Progress: ", i, "/", length(candidates)))
    form <- paste0("cbind(active_count, inactive_count) ~ `", f, "` + ", fixed_f, " + (1 | ", RANDOM_EFF, ")")
    try({
      m <- glmer(as.formula(form), data = df, family = binomial, 
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=1e5)))
      coef <- summary(m)$coefficients
      row_idx <- which(grepl(f, rownames(coef), fixed = TRUE))
      if(length(row_idx) == 0) row_idx <- 2
      
      glmm_res <- rbind(glmm_res, data.frame(
        Feature = f,
        GLMM_Est = coef[row_idx, 1],
        GLMM_P = coef[row_idx, 4],
        Singular = isSingular(m)
      ))
    }, silent = TRUE)
  }
  
  final <- glm_res %>%
    filter(Feature %in% candidates) %>%
    left_join(glmm_res, by = "Feature") %>%
    mutate(GLMM_FDR = p.adjust(GLMM_P, method = "BH")) %>%
    arrange(GLMM_FDR)
    
  return(final)
}

link_annot <- pc_map_linked %>%
  inner_join(annot, by = c("rep_seq_id" = "target")) %>%
  filter(best_hit != "" & !grepl("hypothetical", product, ignore.case = TRUE)) %>%
  select(vOTU, Feature = best_hit, Phage_Name)

mat_annot <- get_prev_matrix(link_annot, otu_stats)
res_annot <- run_two_step_analysis(mat_annot, otu_stats, "BestHit")

if(!is.null(res_annot)) {
  res_annot <- res_annot %>% left_join(besthit_dict, by = "Feature")
  write.csv(res_annot, paste0(OUT_PREFIX, "BestHit.csv"), row.names = FALSE)
}

-
link_pc <- pc_map_linked %>% select(vOTU, Feature = pc_id, Phage_Name)
mat_pc <- get_prev_matrix(link_pc, otu_stats)
res_pc <- run_two_step_analysis(mat_pc, otu_stats, "PC")

if(!is.null(res_pc)) {
  write.csv(res_pc, paste0(OUT_PREFIX, "PC.csv"), row.names = FALSE)
}

