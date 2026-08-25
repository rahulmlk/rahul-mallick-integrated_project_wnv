# Transcriptomics: RNA-seq Differential Expression (DESeq2 + GSVA)

Representative R analysis reflecting the RNA-seq differential expression
workflow used in my thesis for WNV-infected vs. mock-treated macrophages.

## Real pipeline (thesis)

1. Reads aligned to **GRCh38** with **STAR (v2.7)**
2. Gene-level counts via **featureCounts** against **GENCODE v42**
3. Differential expression with **DESeq2 (v1.36)**, Benjamini–Hochberg
   correction, |log2FC| and adjusted p-value thresholds
4. Pathway-level interferon-response scoring with **GSVA**
5. Volcano plot of differentially expressed genes

The upstream alignment/quantification steps (STAR, featureCounts) require
a reference genome, annotation, and raw FASTQ files, so they're
represented here as a documented Bash snippet rather than runnable code.
`example_DEG_analysis.R` starts from a count matrix (as DESeq2 normally
does) and reproduces steps 3-5 on a small synthetic dataset.

## Files

- `example_DEG_analysis.R` - DESeq2 + GSVA + volcano plot, from a count matrix
- `data/sample_counts.csv` - synthetic gene x sample count matrix (20 genes, 8 samples)
- `data/sample_metadata.csv` - synthetic sample condition labels

## Upstream alignment (for reference, not run here)

```bash
# Alignment
STAR --runThreadN 8 --genomeDir GRCh38_index \
     --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix results/sample_

# Quantification
featureCounts -T 8 -p -a gencode.v42.annotation.gtf \
     -o counts.txt results/*_Aligned.sortedByCoord.out.bam
```

## Running the R script

```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GSVA"))
install.packages("ggplot2")
```

```bash
Rscript example_DEG_analysis.R
```

Outputs a results table (`results/deg_results.csv`) and a volcano plot
(`../figures/volcano_plot.png`).
