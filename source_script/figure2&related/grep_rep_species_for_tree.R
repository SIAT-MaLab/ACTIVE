
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
stats_file  <- "prophage_stats_final.tsv"
quality_file <- "allcheckm2_quality.tsv"
genome_dir  <- "/home/huangyan/experiment/SPA/3/gut_isolate/allgenome"
output_dir  <- "species"

df_stats <- read.delim(stats_file, sep = "\t", check.names = FALSE) %>%
  select(Genome_ID, Species) %>%
  filter(!is.na(Species), Species != "")

df_quality <- read.delim(quality_file, sep = "\t", header = FALSE) %>%
  select(Genome_ID = V1, N50 = V7)

representative_df <- df_stats %>%
  inner_join(df_quality, by = "Genome_ID") %>%
  group_by(Species) %>%
  slice_max(order_by = N50, n = 1, with_ties = FALSE) %>%
  ungroup()

write.table(
  representative_df,
  file = "species_best_genomes_by_N50.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

apply(representative_df, 1, function(row) {
  genome_id <- trimws(row["Genome_ID"])
  species   <- trimws(row["Species"])

  species_name <- gsub("\\s+", "_", species)

  src_path  <- file.path(genome_dir, paste0(genome_id, ".fasta"))
  dest_path <- file.path(output_dir, paste0(species_name, ".fasta"))

  if (file.exists(src_path)) {

    if (file.exists(dest_path)) {
      file.remove(dest_path)
    }
    file.symlink(from = normalizePath(src_path), to = dest_path)
  } else {
    warning(paste("missing", src_path))
  }
})


