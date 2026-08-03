#!/usr/bin/env Rscript

# This script runs DESeq2 on our RNA-seq counts data
# We have 4 samples: 2 placebo (pre/post) and 2 mifepristone (pre/post)
# Since we only have 1 sample per group, we cannot let DESeq2 calculate
# dispersion on its own. So we set it manually to 0.4^2 (BCV = 0.4)

# load the DESeq2 library
library(DESeq2)

# get the input arguments (counts file and output folder)
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]
output_folder <- args[2]

# create output folder if it does not exist
dir.create(output_folder, showWarnings = FALSE)

# read the featureCounts output file
# skip = 1 because the first line is just a comment line
counts_data <- read.table(counts_file, header = TRUE, skip = 1, row.names = 1)

# the first 5 columns are gene info (chr, start, end, strand, length)
# so we only keep column 6 onwards, which are the actual sample counts
count_matrix <- counts_data[ , 6:ncol(counts_data)]
count_matrix <- as.matrix(count_matrix)

# rename the columns to just the sample names
colnames(count_matrix) <- c("SRR31630892", "SRR31630893", "SRR31630896", "SRR31630897")

# now we create the sample information table
# column order is: SRR31630892, SRR31630893, SRR31630896, SRR31630897
sample_group <- c("placebo", "placebo", "mifepristone", "mifepristone")
sample_time <- c("pre", "post", "pre", "post")

sample_info <- data.frame(group = sample_group, time = sample_time)
rownames(sample_info) <- colnames(count_matrix)

# make sure group and time are factors
sample_info$group <- factor(sample_info$group, levels = c("placebo", "mifepristone"))
sample_info$time <- factor(sample_info$time, levels = c("pre", "post"))

# build the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = sample_info,
                               design = ~ group * time)

# remove genes with very low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# calculate size factors (needed for normalization)
dds <- estimateSizeFactors(dds)

# since we don't have replicates, we set the dispersion manually instead
# of using estimateDispersions()
dispersions(dds) <- 0.4^2

# run the actual test
dds <- nbinomWaldTest(dds)

# get results for the interaction term
# this tells us if mifepristone caused a bigger change than placebo did
results_table <- results(dds, name = "groupmifepristone.timepost")
results_table <- as.data.frame(results_table)

# add gene names as a column and sort by adjusted p-value
results_table$gene <- rownames(results_table)
results_table <- results_table[order(results_table$padj), ]

# keep only the columns we need
final_results <- data.frame(
  gene = results_table$gene,
  log2FoldChange = results_table$log2FoldChange,
  baseMean = results_table$baseMean,
  pvalue = results_table$pvalue,
  FDR = results_table$padj
)

# save the main result
write.csv(final_results, file = paste0(output_folder, "/deseq2_results.csv"), row.names = FALSE)

# also save the main effects (just group alone and just time alone)
# these are extra, not the main result, but good to have for comparison
drug_effect <- as.data.frame(results(dds, name = "group_mifepristone_vs_placebo"))
time_effect <- as.data.frame(results(dds, name = "time_post_vs_pre"))

write.csv(drug_effect, file = paste0(output_folder, "/deseq2_drug_main_effect.csv"), row.names = TRUE)
write.csv(time_effect, file = paste0(output_folder, "/deseq2_time_main_effect.csv"), row.names = TRUE)

# save normalized counts too, useful for making plots later
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = paste0(output_folder, "/normalized_counts.csv"))

print("Done! Results saved in the output folder.")
