#!/usr/bin/env Rscript

# ============================================================
# DESeq2 differential expression analysis
# Adapted to the uploaded file: data_DEG_R.csv
# Main contrast: Sm_Hf vs n.c
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# ------------------------------------------------------------
# 1. Input / output paths
# ------------------------------------------------------------
input_counts <- "data_DEG_R.csv"
output_dir   <- "results"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read count matrix
# ------------------------------------------------------------
if (!file.exists(input_counts)) {
  stop("Input file not found: ", input_counts)
}

# Keep default check.names=TRUE so duplicate group names become:
# Sm_Hf, Sm_Hf.1, ..., Sm.1, ..., n.c.9
counts <- read.csv(
  input_counts,
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

counts[is.na(counts)] <- 0

counts <- as.matrix(counts)
mode(counts) <- "numeric"

if (anyNA(counts)) {
  stop("NA values are still present after replacement. Please inspect the input CSV.")
}

# Convert to numeric matrix
counts <- as.matrix(counts)
mode(counts) <- "numeric"

if (anyNA(counts)) {
  stop("The count matrix contains NA values. Please inspect the input CSV.")
}

if (nrow(counts) == 0 || ncol(counts) == 0) {
  stop("The count matrix is empty.")
}

# ------------------------------------------------------------
# 3. Derive sample groups from real column names
# ------------------------------------------------------------
sample_names <- colnames(counts)

# Remove suffixes like .1, .2, .3 to recover group labels
condition_raw <- sub("\\.[0-9]+$", "", sample_names)

valid_groups <- c("Sm_Hf", "Sm", "Hf", "n.c")

if (!all(condition_raw %in% valid_groups)) {
  stop(
    "Unexpected group names detected in column names.\n",
    "Detected groups: ", paste(sort(unique(condition_raw)), collapse = ", "), "\n",
    "Expected groups: ", paste(valid_groups, collapse = ", ")
  )
}

condition <- factor(condition_raw, levels = c("n.c", "Hf", "Sm", "Sm_Hf"))

coldata <- data.frame(
  row.names = sample_names,
  condition = condition
)

# Print sample distribution
cat("Sample distribution by group:\n")
print(table(coldata$condition))

# ------------------------------------------------------------
# 4. Filter low-count genes
# ------------------------------------------------------------
# Same rule as in the original script
counts_filtered <- counts[rowSums(counts) > 10, , drop = FALSE]

if (nrow(counts_filtered) == 0) {
  stop("No genes remaining after filtering with rowSums > 10.")
}

# ------------------------------------------------------------
# 5. Build DESeq2 object
# ------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered),
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

# Variance stabilizing transformation for PCA/heatmap
vsd <- vst(dds, blind = FALSE)

# ------------------------------------------------------------
# 6. Differential expression: Sm_Hf vs n.c
# ------------------------------------------------------------
res <- results(dds, contrast = c("condition", "Sm_Hf", "n.c"))
res <- res[order(res$padj), ]

res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_df <- res_df[, c("GeneID", setdiff(colnames(res_df), "GeneID"))]

sig_df <- subset(res_df, !is.na(padj) & padj < 0.05)
sig_lfc_df <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)

# ------------------------------------------------------------
# 7. Save results tables
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 8. PCA plot
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 9. Dispersion plot
# ------------------------------------------------------------
png(
  filename = file.path(output_dir, "dispersion_estimates.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plotDispEsts(dds)
dev.off()

# ------------------------------------------------------------
# 10. MA plot
# ------------------------------------------------------------
png(
  filename = file.path(output_dir, "ma_plot_sm_hf_vs_nc.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plotMA(res, ylim = c(-5, 5), main = "MA plot: Sm_Hf vs n.c")
dev.off()

# ------------------------------------------------------------
# 11. Volcano plot
# ------------------------------------------------------------
volcano_df <- res_df
volcano_df$negLog10P <- -log10(volcano_df$pvalue)
volcano_df$negLog10P[!is.finite(volcano_df$negLog10P)] <- NA

volcano_df$category <- "Not significant"
volcano_df$category[!is.na(volcano_df$padj) & volcano_df$padj < 0.05] <- "padj < 0.05"
volcano_df$category[!is.na(volcano_df$padj) &
                      volcano_df$padj < 0.05 &
                      abs(volcano_df$log2FoldChange) > 1.5] <- "padj < 0.05 & |log2FC| > 1.5"

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10P)) +
  geom_point(aes(color = category), alpha = 0.7, size = 1.8, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey70",
      "padj < 0.05" = "blue",
      "padj < 0.05 & |log2FC| > 1.5" = "red"
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

# ------------------------------------------------------------
# 12. Heatmap
# ------------------------------------------------------------
vsd_mat <- assay(vsd)

if (nrow(sig_df) >= 2) {
  top_genes <- head(sig_df$GeneID, 50)
  top_genes <- top_genes[top_genes %in% rownames(vsd_mat)]
  heatmap_mat <- vsd_mat[top_genes, , drop = FALSE]
} else {
  gene_vars <- apply(vsd_mat, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(50, length(gene_vars))]
  heatmap_mat <- vsd_mat[top_genes, , drop = FALSE]
}

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

# ------------------------------------------------------------
# 13. Save normalized counts
# ------------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$GeneID <- rownames(norm_counts_df)
norm_counts_df <- norm_counts_df[, c("GeneID", setdiff(colnames(norm_counts_df), "GeneID"))]

write.csv(
  norm_counts_df,
  file = file.path(output_dir, "normalized_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Save transformed matrix used in PCA / heatmap
# ------------------------------------------------------------
vsd_df <- as.data.frame(assay(vsd))
vsd_df$GeneID <- rownames(vsd_df)
vsd_df <- vsd_df[, c("GeneID", setdiff(colnames(vsd_df), "GeneID"))]

write.csv(
  vsd_df,
  file = file.path(output_dir, "vst_transformed_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Save session info
# ------------------------------------------------------------
writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

# ------------------------------------------------------------
# 16. Console summary
# ------------------------------------------------------------
cat("\nDESeq2 analysis completed successfully.\n")
cat("Input file: ", input_counts, "\n", sep = "")
cat("Genes before filtering: ", nrow(counts), "\n", sep = "")
cat("Genes after filtering: ", nrow(counts_filtered), "\n", sep = "")
cat("Samples: ", ncol(counts), "\n", sep = "")
cat("Contrast: Sm_Hf vs n.c\n")
cat("Total DE results: ", nrow(res_df), "\n", sep = "")
cat("Significant genes (padj < 0.05): ", nrow(sig_df), "\n", sep = "")
cat("Significant genes (padj < 0.05 and |log2FC| >= 1): ", nrow(sig_lfc_df), "\n", sep = "")
cat("Results written to: ", output_dir, "\n", sep = "")


