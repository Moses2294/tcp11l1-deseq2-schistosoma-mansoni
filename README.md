# tcp11l1-deseq2-schistosoma-mansoni
R/DESeq2 workflow for differential expression analysis of the RNA-seq discovery cohort in school-aged children from rural Cameroon, including the Sm_Hf vs n.c (KK+US+ vs KK-US-) comparison from the TCP11L1 schistosomiasis biomarker study.


# tcp11l1-deseq2-schistosoma-mansoni

This repository contains the R/DESeq2 code used for differential expression analysis of the RNA-seq discovery cohort from the manuscript:

**“All T-complex 11 like 1 is a male-specific blood transcript biomarker of S. mansoni infection in school-aged children”**

## Overview

This repository currently contains the DESeq2-based downstream analysis used in the RNA-seq discovery phase of the study. The broader study investigated host blood transcriptomic changes associated with *Schistosoma mansoni* infection in school-aged children from a rural endemic area in Cameroon.

According to the study design, RNA-seq libraries were generated from four epidemiological groups:
- `Sm_Hf` (KK+US+)
- `Sm` (KK+US-)
- `Hf` (KK-US+)
- `n.c` (KK-US-)

The present script focuses on the comparison:

- **Sm_Hf vs n.c** (KK+US+ vs KK-US-)

## What the script does

The DESeq2 workflow currently:
1. Imports a gene-count matrix (`data_DEG_R.csv`)
2. Filters out very low-count genes (`rowSums > 10`)
3. Creates sample metadata with four conditions
4. Builds a `DESeqDataSet`
5. Runs differential expression analysis with `DESeq2`
6. Applies variance stabilizing transformation (`vst`)
7. Generates exploratory plots:
   - PCA
   - dispersion estimates
   - MA plot
   - volcano-style visualization
8. Extracts the contrast:
   - `Sm_Hf` versus `n.c`
9. Exports:
   - full differential expression results
   - significant genes with adjusted p-value < 0.05
## Main script

The main analysis script is:

```text
#!/usr/bin/env Rscript

# ============================================================
# DESeq2 differential expression analysis
# Contrast: Sm_Hf vs n.c
# Study context: Schistosoma mansoni RNA-seq discovery cohort
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# -----------------------------
# 1. User-defined input/output
# -----------------------------
input_counts <- "data_DEG_R.csv"
output_dir   <- "results"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Read count matrix
# -----------------------------
if (!file.exists(input_counts)) {
  stop("Input file not found: ", input_counts)
}

counts <- read.csv(
  input_counts,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Convert to numeric matrix
counts <- as.matrix(counts)
mode(counts) <- "numeric"

if (any(is.na(counts))) {
  stop("Count matrix contains NA values. Please check the input file.")
}

if (ncol(counts) == 0 || nrow(counts) == 0) {
  stop("Count matrix is empty.")
}

# ------------------------------------------
# 3. Define sample groups (original ordering)
# ------------------------------------------
# Original script assumes 40 samples ordered as:
# 10 Sm_Hf, 10 Sm, 10 Hf, 10 n.c

expected_n <- 40
if (ncol(counts) != expected_n) {
  stop(
    "This script expects ", expected_n,
    " columns/samples in the following order:\n",
    "10 Sm_Hf, 10 Sm, 10 Hf, 10 n.c\n",
    "Detected columns: ", ncol(counts), "\n",
    "Please adapt the condition vector if your sample order differs."
  )
}

condition <- factor(
  c(
    rep("Sm_Hf", 10),
    rep("Sm",    10),
    rep("Hf",    10),
    rep("n.c",   10)
  ),
  levels = c("n.c", "Hf", "Sm", "Sm_Hf")
)

coldata <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

# -------------------------------------
# 4. Filter low-count genes
# -------------------------------------
# Faithful to original script: keep genes with rowSums > 10
counts_filtered <- counts[rowSums(counts) > 10, ]

if (nrow(counts_filtered) == 0) {
  stop("No genes remaining after filtering (rowSums > 10).")
}

# -------------------------------------
# 5. Build DESeq2 dataset and run model
# -------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

# Variance stabilizing transformation for exploratory plots
vsd <- vst(dds, blind = FALSE)

# -------------------------------------
# 6. Main contrast: Sm_Hf vs n.c
# -------------------------------------
res <- results(dds, contrast = c("condition", "Sm_Hf", "n.c"))
res <- res[order(res$padj), ]

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

sig_df <- subset(res_df, !is.na(padj) & padj < 0.05)

# Optional stricter subset with effect size threshold
sig_lfc_df <- subset(sig_df, abs(log2FoldChange) >= 1)

# -------------------------------------
# 7. Save result tables
# -------------------------------------
write.csv(
  res_df,
  file = file.path(output_dir, "deseq2_full_results_sm_hf_vs_nc.csv"),
  row.names = FALSE
)

write.csv(
  sig_df,
  file = file.path(output_dir, "deseq2_significant_results_padj_lt_0.05_sm_hf_vs_nc.csv"),
  row.names = FALSE
)

write.csv(
  sig_lfc_df,
  file = file.path(output_dir, "deseq2_significant_results_padj_lt_0.05_absLFC_ge_1_sm_hf_vs_nc.csv"),
  row.names = FALSE
)

# -------------------------------------
# 8. PCA plot
# -------------------------------------
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of variance-stabilized counts") +
  theme_bw()

ggsave(
  filename = file.path(output_dir, "pca_sm_hf_vs_nc.png"),
  plot = p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------------------
# 9. Dispersion plot
# -------------------------------------
png(
  filename = file.path(output_dir, "dispersion_estimates.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plotDispEsts(dds)
dev.off()

# -------------------------------------
# 10. MA plot
# -------------------------------------
png(
  filename = file.path(output_dir, "ma_plot_sm_hf_vs_nc.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plotMA(res, ylim = c(-5, 5), main = "MA plot: Sm_Hf vs n.c")
dev.off()

# -------------------------------------
# 11. Volcano plot
# -------------------------------------
volcano_df <- res_df
volcano_df$negLog10P <- -log10(volcano_df$pvalue)
volcano_df$category <- "Not significant"
volcano_df$category[!is.na(volcano_df$padj) & volcano_df$padj < 0.01] <- "padj < 0.01"
volcano_df$category[!is.na(volcano_df$padj) &
                      volcano_df$padj < 0.01 &
                      abs(volcano_df$log2FoldChange) > 2] <- "padj < 0.01 & |log2FC| > 2"

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10P)) +
  geom_point(aes(color = category), alpha = 0.7, size = 1.8) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey70",
      "padj < 0.01" = "blue",
      "padj < 0.01 & |log2FC| > 2" = "red"
    )
  ) +
  labs(
    title = "Volcano plot: Sm_Hf vs n.c",
    x = "log2 fold change",
    y = "-log10(p-value)",
    color = NULL
  ) +
  theme_bw()

ggsave(
  filename = file.path(output_dir, "volcano_sm_hf_vs_nc.png"),
  plot = p_volcano,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------------------
# 12. Heatmap
# -------------------------------------
# Use top significant genes if available; otherwise top variable genes
vsd_mat <- assay(vsd)

if (nrow(sig_df) >= 2) {
  top_genes <- head(sig_df$gene, 50)
  top_genes <- top_genes[top_genes %in% rownames(vsd_mat)]
  heatmap_mat <- vsd_mat[top_genes, , drop = FALSE]
} else {
  gene_vars <- apply(vsd_mat, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(50, length(gene_vars))]
  heatmap_mat <- vsd_mat[top_genes, , drop = FALSE]
}

# Row scaling for visualization
heatmap_mat <- t(scale(t(heatmap_mat)))
heatmap_mat[is.na(heatmap_mat)] <- 0

annotation_col <- data.frame(condition = coldata$condition)
rownames(annotation_col) <- rownames(coldata)

png(
  filename = file.path(output_dir, "heatmap_top_genes_sm_hf_vs_nc.png"),
  width = 2000,
  height = 1800,
  res = 220
)
pheatmap(
  heatmap_mat,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  main = "Heatmap of top genes"
)
dev.off()

# -------------------------------------
# 13. Save normalized counts
# -------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene <- rownames(norm_counts_df)
norm_counts_df <- norm_counts_df[, c("gene", setdiff(colnames(norm_counts_df), "gene"))]

write.csv(
  norm_counts_df,
  file = file.path(output_dir, "normalized_counts.csv"),
  row.names = FALSE
)

# -------------------------------------
# 14. Save session information
# -------------------------------------
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

# -------------------------------------
# 15. Console summary
# -------------------------------------
cat("DESeq2 analysis completed successfully.\n")
cat("Input file: ", input_counts, "\n", sep = "")
cat("Genes before filtering: ", nrow(counts), "\n", sep = "")
cat("Genes after filtering: ", nrow(counts_filtered), "\n", sep = "")
cat("Total DE results: ", nrow(res_df), "\n", sep = "")
cat("Significant genes (padj < 0.05): ", nrow(sig_df), "\n", sep = "")
cat("Significant genes (padj < 0.05 and |log2FC| >= 1): ", nrow(sig_lfc_df), "\n", sep = "")
cat("Results written to: ", output_dir, "\n", sep = "")
