#!/usr/bin/env Rscript

library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
output_folder <- args[2]

dir.create(output_folder, showWarnings = FALSE)

counts_data <- read.table(counts_file, header = TRUE, skip = 1, row.names = 1)

count_matrix <- counts_data[ , 6:ncol(counts_data)]
count_matrix <- as.matrix(count_matrix)

colnames(count_matrix) <- c("SRR31630892", "SRR31630893", "SRR31630896", "SRR31630897")

sample_group <- c("placebo", "placebo", "mifepristone", "mifepristone")
sample_time <- c("pre", "post", "pre", "post")

sample_info <- data.frame(group = sample_group, time = sample_time)
rownames(sample_info) <- colnames(count_matrix)

sample_info$group <- factor(sample_info$group, levels = c("placebo", "mifepristone"))
sample_info$time <- factor(sample_info$time, levels = c("pre", "post"))

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = sample_info,
                               design = ~ group * time)

# remove genes with very low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds <- estimateSizeFactors(dds)


dispersions(dds) <- 0.4^2

dds <- nbinomWaldTest(dds)

results_table <- results(dds, name = "groupmifepristone.timepost")
results_table <- as.data.frame(results_table)

results_table$gene <- rownames(results_table)
results_table <- results_table[order(results_table$padj), ]

final_results <- data.frame(
  gene = results_table$gene,
  log2FoldChange = results_table$log2FoldChange,
  baseMean = results_table$baseMean,
  pvalue = results_table$pvalue,
  FDR = results_table$padj
)

write.csv(final_results, file = paste0(output_folder, "/deseq2_results.csv"), row.names = FALSE)

drug_effect <- as.data.frame(results(dds, name = "group_mifepristone_vs_placebo"))
time_effect <- as.data.frame(results(dds, name = "time_post_vs_pre"))

write.csv(drug_effect, file = paste0(output_folder, "/deseq2_drug_main_effect.csv"), row.names = TRUE)
write.csv(time_effect, file = paste0(output_folder, "/deseq2_time_main_effect.csv"), row.names = TRUE)

normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = paste0(output_folder, "/normalized_counts.csv"))

print("Done! Results saved in the output folder.")
