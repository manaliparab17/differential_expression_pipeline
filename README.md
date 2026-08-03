RNA-seq Differential Expression: Mifepristone vs Placebo (BRCA carriers)

RNA-seq pipeline for a 2x2 design (drug/placebo x pre/post) comparing endometrial gene expression before and after short-term mifepristone treatment in BRCA1/2 mutation carriers. Raw reads -> QC -> alignment -> gene counts -> DESeq2.

Note: there is only 1 sample per group (no replicates). This means DESeq2 can't estimate dispersion normally, so I set it manually (BCV = 0.4) in deseq2_analysis.R. The main test used is the group:time interaction term, which basically asks "did expression change more after mifepristone than after placebo?" (the placebo arm is there to control for changes that happen just from time/handling, not the drug itself).
What the pipeline does
SRA reads -> FastQC -> fastp -> HISAT2 -> samtools sort/index
-> featureCounts -> DESeq2

How to set it up
bash
conda create -n diffexp -c conda-forge -c bioconda -c defaults \
sra-tools \
fastqc \
fastp \
hisat2 \
samtools \
subread \
"r-base>=4.0" \
bioconductor-deseq2

conda activate diffexp

Building the HISAT2 index from the full genome needs a lot of RAM (~160GB), so if your machine can't handle that, you can download an already-built index instead:

bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
How to run it
bash
chmod +x differential_expression.sh
./differential_expression.sh

I also added two scripts to check the results visually before trusting them (this actually caught a real problem in my data):

bash
Rscript pca_plot.R normalized_counts.csv DE/
Rscript diagnostic_plots.R deseq2_results.csv normalized_counts.csv DE/
