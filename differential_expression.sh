#!/bin/bash

set -e

mkdir -p Sequences Reference QC Trimmed Alignment Counts DE

# download reads
cd Sequences

for SAMPLE in SRR31630892 SRR31630893 SRR31630896 SRR31630897; do
    prefetch $SAMPLE
    fasterq-dump --split-files $SAMPLE -O .
    gzip ${SAMPLE}_*.fastq
done

cd ..

SAMPLES=(SRR31630892 SRR31630893 SRR31630896 SRR31630897)

# reference genome + annotation (GRCh38)

cd Reference

wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa

wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz
mv Homo_sapiens.GRCh38.110.gtf annotation.gtf

cd ..

hisat2-build Reference/genome.fa Reference/grch38_index
samtools faidx Reference/genome.fa

fastqc Sequences/*.fastq.gz -o QC

for SAMPLE in "${SAMPLES[@]}"; do
    fastp -i Sequences/${SAMPLE}_1.fastq.gz -I Sequences/${SAMPLE}_2.fastq.gz \
        -o Trimmed/${SAMPLE}_1.trimmed.fastq.gz -O Trimmed/${SAMPLE}_2.trimmed.fastq.gz \
        -j QC/${SAMPLE}_fastp.json -h QC/${SAMPLE}_fastp.html
done

for SAMPLE in "${SAMPLES[@]}"; do
    hisat2 -p 4 -x Reference/grch38_index \
        -1 Trimmed/${SAMPLE}_1.trimmed.fastq.gz \
        -2 Trimmed/${SAMPLE}_2.trimmed.fastq.gz \
        -S Alignment/${SAMPLE}.sam
done

for SAMPLE in "${SAMPLES[@]}"; do
    samtools view -Sb Alignment/${SAMPLE}.sam > Alignment/${SAMPLE}.bam
    samtools sort Alignment/${SAMPLE}.bam -o Alignment/${SAMPLE}.sorted.bam
    samtools index Alignment/${SAMPLE}.sorted.bam
    samtools flagstat Alignment/${SAMPLE}.sorted.bam > Alignment/${SAMPLE}_alignment_stats.txt
done

featureCounts -T 4 -p --countReadPairs -a Reference/annotation.gtf \
    -o Counts/gene_counts.txt \
    Alignment/SRR31630892.sorted.bam \
    Alignment/SRR31630893.sorted.bam \
    Alignment/SRR31630896.sorted.bam \
    Alignment/SRR31630897.sorted.bam

Rscript deseq2_analysis.R Counts/gene_counts.txt DE/

# significant DEGs: FDR < 0.05, |log2FC| > 1
awk -F',' 'NR==1 || ($5 != "NA" && $5+0 < 0.05 && ($2+0 > 1 || $2+0 < -1))' \
    DE/deseq2_results.csv > DE/significant_DEGs.csv

echo "done"
