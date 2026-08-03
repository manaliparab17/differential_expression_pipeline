#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
output_folder <- args[2]

dir.create(output_folder, showWarnings = FALSE)

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)

sample_names <- colnames(counts)
sample_labels <- c("placebo_pre", "placebo_post", "mife_pre", "mife_post")
sample_colors <- c("#2a78d6", "#2a78d6", "#eb6834", "#eb6834")

log_counts <- log2(counts + 1)


pca <- prcomp(t(log_counts))

var_explained <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 1)

png(file.path(output_folder, "pca_plot.png"), width = 700, height = 600)

plot(pca$x[, 1], pca$x[, 2],
     col = sample_colors, pch = 19, cex = 3,
     xlab = paste0("PC1 (", var_explained[1], "% variance)"),
     ylab = paste0("PC2 (", var_explained[2], "% variance)"),
     main = "PCA of samples")

text(pca$x[, 1], pca$x[, 2], labels = sample_labels, pos = 3, cex = 0.9)

legend("topright", legend = c("placebo", "mifepristone"),
       col = c("#2a78d6", "#eb6834"), pch = 19)

dev.off()

lib_sizes <- colSums(counts)

png(file.path(output_folder, "library_size.png"), width = 700, height = 500)

barplot(lib_sizes,
        names.arg = sample_labels,
        col = sample_colors,
        main = "Total normalized library size per sample",
        ylab = "Total normalized counts",
        las = 2)

dev.off()

print("Saved pca_plot.png and library_size.png")
