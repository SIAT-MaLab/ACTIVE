suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})


.parse_dir <- function(default = ".") {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep("^--dir=", args, value = TRUE)
  if (length(hit) > 0) return(sub("^--dir=", "", hit[1]))
  idx <- which(args == "--dir")
  if (length(idx) > 0 && length(args) >= idx + 1) return(args[idx + 1])
  return(default)
}


target_dir <- .parse_dir()
csv_files <- list.files(path = target_dir, pattern = "\\.csv$", full.names = TRUE)


for (f in csv_files) {

  df <- tryCatch(read.csv(f), error = function(e) NULL)

  
  deltas <- df$Cliffs_Delta
  n_samples <- length(deltas)
  
  if (n_samples < 3) {
    cat(sprintf("n_samples < 3\n", basename(f), n_samples))
    next
  }
  

  base_name <- tools::file_path_sans_ext(basename(f))
  match <- str_match(base_name, "^(.*)_(\\d+)$")
  if (!is.na(match[1,1])) {
    phage_name <- match[1,2]
    n_iters <- match[1,3]
  } else {
    phage_name <- base_name
    n_iters <- as.character(n_samples)
  }


  mean_val <- mean(deltas, na.rm = TRUE)
  ci_lower <- quantile(deltas, 0.025, na.rm = TRUE, names = FALSE)
  ci_upper <- quantile(deltas, 0.975, na.rm = TRUE, names = FALSE)
  sw_test <- shapiro.test(deltas)
  p_val <- sw_test$p.value
  
  p_text <- ifelse(p_val < 0.0001, "< 0.0001", sprintf("%.4f", p_val))
  subtitle_text <- sprintf("Mean = %.4f | 95%% CI = [%.4f, %.4f]\nShapiro-Wilk p-value: %s", 
                           mean_val, ci_lower, ci_upper, p_text)
  p <- ggplot(df, aes(x = Cliffs_Delta)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, 
                   fill = "#84AADD", color = "white", alpha = 0.9) +
    geom_density(color = "#D94250", linewidth = 1) +
    geom_vline(xintercept = mean_val, color = "#1D2883", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = ci_lower, color = "#F6A834", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = ci_upper, color = "#F6A834", linetype = "dotted", linewidth = 1) +
    labs(
      title = sprintf("Cliff's Delta Sampling Distribution (%s, n = %s)", phage_name, n_iters),
      subtitle = subtitle_text,
      x = "Cliff's Delta",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#EBEBEB"),
      panel.grid.major.y = element_line(color = "#EBEBEB"),
      axis.line.x = element_line(color = "black")
    )
  
  out_pdf <- file.path(target_dir, sprintf("%s_distribution.pdf", base_name))
  ggsave(filename = out_pdf, plot = p, width = 7, height = 5, bg = "white")
  
  cat(sprintf("  -> %s (p = %s)\n", basename(out_pdf), p_text))
}

