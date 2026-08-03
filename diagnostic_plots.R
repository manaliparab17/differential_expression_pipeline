#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
results_file <- args[1]
counts_file <- args[2]
output_folder <- args[3]

dir.create(output_folder, showWarnings = FALSE)

results <- read.csv(results_file)

fdr_cutoff <- 0.05
fc_cutoff <- 1

results$is_significant <- !is.na(results$FDR) &
  results$FDR < fdr_cutoff &
  abs(results$log2FoldChange) > fc_cutoff

results$neg_log10_fdr <- -log10(results$FDR)

point_colors <- ifelse(results$is_significant, "#e34948", "#b4b2a9")

png(file.path(output_folder, "volcano_plot.png"), width = 700, height = 600)

plot(results$log2FoldChange, results$neg_log10_fdr,
     col = point_colors, pch = 19, cex = 0.6,
     xlab = "log2 fold change",
     ylab = "-log10(FDR)",
     main = "Volcano plot")

abline(v = c(-fc_cutoff, fc_cutoff), lty = 2, col = "gray40")
abline(h = -log10(fdr_cutoff), lty = 2, col = "gray40")

legend("topright", legend = c("significant", "not significant"),
       col = c("#e34948", "#b4b2a9"), pch = 19)

dev.off()

results$log_basemean <- log2(results$baseMean + 1)

png(file.path(output_folder, "ma_plot.png"), width = 700, height = 600)

plot(results$log_basemean, results$log2FoldChange,
     col = point_colors, pch = 19, cex = 0.6,
     xlab = "log2 mean expression",
     ylab = "log2 fold change",
     main = "MA plot")

abline(h = 0, lty = 2, col = "gray40")

legend("topright", legend = c("significant", "not significant"),
       col = c("#e34948", "#b4b2a9"), pch = 19)

dev.off()

norm_counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
norm_counts <- as.matrix(norm_counts)

sig_genes <- results[results$is_significant, ]
sig_genes <- sig_genes[order(sig_genes$FDR), ]
top_n <- min(20, nrow(sig_genes))
top_genes <- sig_genes$gene[1:top_n]

if (top_n > 1) {

  heatmap_data <- norm_counts[top_genes, ]

  heatmap_data <- t(scale(t(log2(heatmap_data + 1))))

  colnames(heatmap_data) <- c("placebo_pre", "placebo_post", "mife_pre", "mife_post")

  png(file.path(output_folder, "heatmap.png"), width = 700, height = 800)

  heatmap(heatmap_data,
          scale = "none",
          col = colorRampPalette(c("#378ADD", "white", "#E24B4A"))(50),
          margins = c(8, 10),
          main = "Top significant genes")

  dev.off()

} else {
  print("Not enough significant genes to make a heatmap")
}

print("Saved volcano_plot.png, ma_plot.png, and heatmap.png (if enough genes)")
