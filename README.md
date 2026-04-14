# tcp11l1-deseq2-schistosoma-mansoni
R/DESeq2 workflow for differential expression analysis of the RNA-seq discovery cohort in school-aged children from rural Cameroon, including the Sm_Hf vs n.c comparison from the TCP11L1 schistosomiasis biomarker study.


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

- **Sm_Hf vs n.c**

## What the script does

The DESeq2 workflow currently:
1. Imports a gene-count matrix (`DEG_R.csv`)
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

## Repository structure

A recommended structure is:

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── scripts/
│   └── deseq2_sm_hf_vs_nc.R
├── data/
│   └── DEG_R.csv
└── results/
    ├── deseq2_full_results_sm_hf_vs_nc.csv
    └── deseq2_significant_results_padj_0.05_sm_hf_vs_nc.csv
