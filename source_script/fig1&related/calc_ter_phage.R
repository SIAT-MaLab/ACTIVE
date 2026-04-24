suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})


gbic_fasta     <- "/home/huangyan/experiment/SPA/2-2/GBIC/ref/gbic.fasta"
gbic_coord     <- "/home/huangyan/experiment/SPA/2-2/GBIC/gbic_coordinate.tsv"
gbic_propagate <- "/home/huangyan/experiment/SPA/2-2/used_results/GBIC_propagate.tsv"
ab_fasta       <- "/home/huangyan/experiment/SPA/2-2/AB/CM_results/genome/CM_AB.fasta"
ab_coord       <- "/home/huangyan/experiment/SPA/2-2/AB/CM_results/CM_coordinate.tsv"
ab_propagate   <- "/home/huangyan/experiment/SPA/2-2/used_results/Con_AB_propagate.tsv"

out_plot       <- "./plot_output/Terminal_Phage_PieChart.pdf"

min_length <- 100  
label_internal    <- "Internal Phage"
label_term_comp   <- "Computable Terminal Phage"
label_term_incomp <- "Incomputable Terminal Phage"
color_internal    <- "grey60"  
color_term_comp   <- "#BF352A" 
color_term_incomp <- "#F6A834" 
pie_border_color <- "black"
pie_border_size  <- 0.7


plot_width       <- 14  
plot_height      <- 6

legend_position  <- "bottom"
legend_text_size <- 13
title_size       <- 16
title_face       <- "bold"

get_fasta_lengths <- function(fasta_file) {
  awk_cmd <- sprintf("awk '/^>/{if(N) print N\"\\t\"L; N=substr($1,2); L=0; next} {L+=length($0)} END{print N\"\\t\"L}' %s", fasta_file)
  dt_lengths <- fread(cmd = awk_cmd, header = FALSE, col.names = c("contig", "contig_length"))
  return(dt_lengths)
}

if (!dir.exists(dirname(out_plot))) dir.create(dirname(out_plot), recursive = TRUE)

process_and_plot_pie <- function(coord_file, fasta_file, propagate_file, source_name) {
  dt_coord <- fread(coord_file, header = FALSE, col.names = c("contig", "start", "stop"))
  dt_lens <- get_fasta_lengths(fasta_file)
  dt_prop <- fread(propagate_file, header = TRUE, select = c("prophage", "host_len"))
  dt_merged <- merge(dt_coord, dt_lens, by = "contig", all.x = TRUE)

  if (any(is.na(dt_merged$contig_length))) {
    warning("  -> warning ", sum(is.na(dt_merged$contig_length)), " loss")
    dt_merged <- dt_merged[!is.na(contig_length)]
  }

  dt_merged[, prophage := paste(contig, start, stop, sep = "_")]
  dt_merged <- merge(dt_merged, dt_prop, by = "prophage", all.x = TRUE)
  dt_merged[, Type := label_internal]

  is_terminal <- dt_merged$start < min_length | (dt_merged$contig_length - dt_merged$stop) < min_length

  dt_merged[is_terminal, Type := ifelse(!is.na(host_len) & host_len == 0, 
                                        label_term_incomp, 
                                        label_term_comp)]


  dt_summary <- dt_merged[, .N, by = Type]
  dt_summary[, fraction := N / sum(N)]
  dt_summary[, label := sprintf("%s\n(n=%d, %.1f%%)", Type, N, fraction * 100)]

  all_types <- data.table(Type = c(label_term_incomp, label_term_comp, label_internal))
  dt_summary <- merge(all_types, dt_summary, by = "Type", all.x = TRUE)
  dt_summary[is.na(N), `:=`(N = 0, fraction = 0, label = "")]

  dt_summary$Type <- factor(dt_summary$Type, levels = c(label_internal, label_term_comp, label_term_incomp))

  p <- ggplot(dt_summary[N > 0], aes(x = "", y = fraction, fill = Type)) +
    geom_bar(stat = "identity", width = 1, color = pie_border_color, linewidth = pie_border_size) +
    coord_polar("y", start = 0) +
    scale_fill_manual(
      values = setNames(c(color_internal, color_term_comp, color_term_incomp), 
                        c(label_internal, label_term_comp, label_term_incomp)),
      labels = setNames(dt_summary$label, dt_summary$Type) 
    ) +
    labs(title = paste0("Phage Location Profile: ", source_name)) +
    theme_void() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = title_face, size = title_size, margin = margin(b = 20)),
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(1.5, "lines")
    )

  return(p)
}


plot_gbic <- process_and_plot_pie(gbic_coord, gbic_fasta, gbic_propagate, "GBIC")
plot_ab   <- process_and_plot_pie(ab_coord, ab_fasta, ab_propagate, "AB")


combined_plot <- (plot_gbic | plot_ab) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = legend_position)


ggsave(out_plot, combined_plot, width = plot_width, height = plot_height, bg = "white")

