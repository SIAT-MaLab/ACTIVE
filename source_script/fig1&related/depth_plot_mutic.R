library(dplyr)
library(ggplot2)

phage_list <- readLines("c.txt")  # phage_list.txt

for (phage_name in phage_list) {
  parts <- strsplit(phage_name, "_")[[1]]
  contig_id <- paste(parts[1:2], collapse = "_")
  start_pos <- as.integer(parts[3])
  end_pos   <- as.integer(parts[4])
  genome_name <- parts[1]

  depth_file_spi <- paste0(genome_name, "_host_depth_spi.txt")
  depth_file_mmc <- paste0(genome_name, "_host_depth_mmc.txt")
  depth_spi <- read.table(depth_file_spi, header = FALSE, sep = "\t", col.names = c("contig", "position", "depth"))
  depth_mmc <- read.table(depth_file_mmc, header = FALSE, sep = "\t", col.names = c("contig", "position", "depth"))
  median_spi <- median(depth_spi$depth, na.rm = TRUE)
  median_mmc <- median(depth_mmc$depth, na.rm = TRUE)
  depth_spi <- depth_spi %>%
    mutate(depth_spi = depth / median_spi)  
  depth_mmc <- depth_mmc %>%
    mutate(depth_mmc = depth / median_mmc)  
  depth_spi$condition <- "spi"
  depth_mmc$condition <- "mmc"
  depth_combined <- full_join(depth_spi, depth_mmc, by = c("contig", "position"))
  depth_combined <- depth_combined %>%
    mutate(contig_number = as.numeric(sub(".*_(\\d+)$", "\\1", contig))) %>%
    arrange(contig_number, position)  
  depth_combined <- depth_combined %>% mutate(global_row = row_number())
  start_row <- depth_combined %>%
    filter(contig == contig_id & position == start_pos) %>%
    pull(global_row) %>%
    head(1) 
  end_row <- depth_combined %>%
    filter(contig == contig_id & position == end_pos) %>%
    pull(global_row) %>%
    tail(1) 

  if (is.na(start_row) | is.na(end_row)) {
    cat("Error: start_row or end_row is NA. Skipping phage: ", phage_name, "\n")
    next  
  }

  if (start_row > end_row) {
    tmp <- start_row
    start_row <- end_row
    end_row <- tmp
  }

  pad <- 30000L
  local_start_row <- max(start_row - pad, 1L)  
  local_end_row   <- min(end_row + pad, nrow(depth_combined))  

  local_region_by_row <- depth_combined[local_start_row:local_end_row, ]

  p_local <- ggplot() +
  geom_line(data = local_region_by_row, aes(x = global_row, y = depth_spi, color = "spi"), alpha = 0.8) +
  geom_line(data = local_region_by_row, aes(x = global_row, y = depth_mmc, color = "mmc"), alpha = 0.8) +
  geom_rect(aes(xmin = 114048, xmax = 128198, ymin = 0, ymax = Inf), 
            fill = "orange", alpha = 0.15) +  
  geom_rect(aes(xmin = 128165, xmax = 210193, ymin = 0, ymax = Inf),
            fill = "gray", alpha = 0.3) +  
  labs(title = paste0(phage_name, " (by global index)"),
       x = "Position (bp)", y = "Normalized Depth") +
  
  theme_bw() +
  scale_color_manual(values = c("spi" = "blue", "mmc" = "red")) 
  ggsave(paste0(phage_name, "_combined.png"), p_local, width = 7, height = 4, dpi = 2400)

  cat("Saved plot for phage: ", phage_name, "\n")  
}

