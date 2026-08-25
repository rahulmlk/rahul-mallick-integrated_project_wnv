#!/usr/bin/env Rscript
#
# example_DEG_analysis.R
#
# Representative differential expression + pathway scoring workflow,
# reflecting the DESeq2/GSVA analysis used in my thesis on WNV-infected
# vs. mock-treated macrophages (STAR/featureCounts -> DESeq2 -> GSVA ->
# volcano plot). Runs on a small synthetic count matrix included in
# data/ so it can be executed standalone without raw sequencing data.
#
# Usage:
#   Rscript example_DEG_analysis.R
#
# Requires: DESeq2, GSVA, ggplot2 (Bioconductor + CRAN)

suppressMessages({
  library(DESeq2)
  library(GSVA)
  library(ggplot2)
})

counts <- read.csv("data/sample_counts.csv", row.names = 1, check.names = FALSE)
meta   <- read.csv("data/sample_metadata.csv", row.names = 1)
meta$condition <- factor(meta$condition, levels = c("mock", "WNV_infected"))

stopifnot(all(colnames(counts) == rownames(meta)))

## ---- 1. Differential expression (DESeq2) --------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~condition
)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "WNV_infected", "mock"),
                alpha = 0.05)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Benjamini-Hochberg correction is already applied by results() (padj);
# flag significance the same way as the thesis pipeline: |log2FC| > 1 and
# padj < 0.05
res_df$significant <- with(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

dir.create("results", showWarnings = FALSE)
write.csv(res_df[order(res_df$padj), ], "results/deg_results.csv", row.names = FALSE)

cat(sprintf("DESeq2: %d of %d genes called significant (padj < 0.05, |log2FC| > 1)\n",
            sum(res_df$significant), nrow(res_df)))

## ---- 2. Pathway-level scoring (GSVA) -------------------------------------
## Mirrors the thesis's interferon-stimulated-gene (ISG) pathway scoring:
## rather than looking at individual genes, score each sample for
## enrichment of a defined gene set.

isg_gene_set <- list(interferon_stimulated_genes = rownames(counts)[1:6])

norm_counts <- counts(dds, normalized = TRUE)
gsva_par <- gsvaParam(as.matrix(log2(norm_counts + 1)), isg_gene_set)
gsva_scores <- gsva(gsva_par)

gsva_df <- data.frame(
  sample = colnames(gsva_scores),
  isg_score = as.numeric(gsva_scores["interferon_stimulated_genes", ]),
  condition = meta[colnames(gsva_scores), "condition"]
)
write.csv(gsva_df, "results/gsva_scores.csv", row.names = FALSE)

cat("\nGSVA interferon-stimulated-gene pathway scores by condition:\n")
print(aggregate(isg_score ~ condition, gsva_df, mean))

## ---- 3. Volcano plot ------------------------------------------------------

res_df$neg_log10_padj <- -log10(res_df$padj)
res_df$neg_log10_padj[is.infinite(res_df$neg_log10_padj)] <- NA

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj,
                               color = significant)) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(values = c("grey70", "firebrick"),
                      labels = c("Not significant", "Significant"),
                      name = NULL) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(title = "Differentially expressed genes: WNV-infected vs. mock",
       x = "log2 fold change", y = "-log10 adjusted p-value") +
  theme_minimal(base_size = 12)

dir.create("../figures", showWarnings = FALSE)
ggsave("../figures/volcano_plot.png", volcano, width = 6, height = 5, dpi = 150)

cat("\nWrote results/deg_results.csv, results/gsva_scores.csv, ../figures/volcano_plot.png\n")
