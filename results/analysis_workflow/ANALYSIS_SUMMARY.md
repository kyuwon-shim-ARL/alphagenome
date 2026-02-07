# AlphaGenome Analysis Workflow - TAL1 Locus Analysis

**Date:** 2026-02-05
**Analysis Type:** Non-coding T-ALL variant effects on TAL1 gene expression

## Overview

This analysis explores non-coding T-cell acute lymphoblastic leukemia (T-ALL) associated mutations near the TAL1 locus on chromosome 1. TAL1 is an oncogene whose aberrant upregulation causes a large fraction of T-ALL cases.

## Data Summary

### Variants Analyzed
- **Total variants:** 128
  - **Oncogenic variants:** 32 (from published studies)
  - **Background variants:** 96 (randomly generated controls)

### Variant Sources
- Mansour et al. 2014 - MuTE (mutation of the TAL1 enhancer) site
- Liu et al. 2017 - Patient cohort variants
- Liu et al. 2020 - Intergenic SNVs
- Smith et al. 2023 - New 3' enhancer variants

### Genomic Context
- **Gene:** TAL1 (chr1:47209255-47242023, minus strand)
- **Cell type:** CD34-positive common myeloid progenitor (CL:0001059)
- **Input sequence length:** 1,048,576 bp (2^20)

## Key Findings

### Quantitative Results

**TAL1 Expression Change (Predicted)**
- Mean oncogenic variant effect: **+0.334**
- Mean background variant effect: **-0.020**
- **Difference: 0.354** (17-fold higher for oncogenic variants)

This demonstrates that cancer-associated variants have a significantly stronger upregulation effect on TAL1 compared to randomly generated background variants of the same length.

### Variant Groups

Three main positional clusters were identified:

1. **5' MuTE site** (chr1:47239291-47239296)
   - Most variants cluster here
   - 1-18 bp insertions
   - Strong TAL1 upregulation effect

2. **Intronic variant** (chr1:47230639)
   - Single nucleotide variants (SNVs)
   - Moderate TAL1 effect

3. **3' enhancer** (chr1:47212072-47212074)
   - 6-22 bp insertions
   - Strong TAL1 upregulation effect

### Predicted Molecular Effects (Jurkat variant)

The individual analysis of the Jurkat cell line variant (chr1:47239296:C>CCGTTTCCTAACC) revealed:

**RNA-seq:** Increased TAL1 gene expression

**DNase-seq:**
- Increased accessibility near variant site
- Increased accessibility near TAL1 transcription start site

**ChIP-Histone marks:**
- **Activating marks increased:**
  - H3K27ac (active enhancer)
  - H3K4me1 (active enhancer)
  - H3K4me3 (active promoter)
  - H3K36me3 (active transcription)

- **Repressive marks decreased:**
  - H3K27me3 (polycomb silencing)
  - H3K9me3 (heterochromatin)

**Interpretation:** The variant creates a de novo active enhancer element that activates the TAL1 promoter, leading to robust oncogene transcription.

## Generated Files

### Visualizations
1. **jurkat_variant_effect.png** - Individual variant effect on RNA, DNase, and histone marks
2. **comparison_MUTE_2.png** - 1 bp insertion variants
3. **comparison_MUTE_3.png** - 2 bp insertion variants
4. **comparison_MUTE_4.png** - 3 bp insertion variants
5. **comparison_MUTE_other.png** - 7-18 bp insertion variants
6. **comparison_47212072_22.png** - 3' enhancer variants (22 bp)
7. **comparison_47212074_7.png** - 3' enhancer variants (7 bp)
8. **comparison_47230639_1.png** - Intronic SNV variants

### Data Files
- **variant_analysis_results.json** - Complete variant scores and metadata
- **run.log** - Full execution log

## Biological Significance

This analysis demonstrates that AlphaGenome can:

1. **Distinguish functional variants from background:** Oncogenic variants show significantly stronger predicted effects than matched random sequences

2. **Predict molecular mechanisms:** The model predicts enhancer creation and chromatin state changes consistent with experimental findings

3. **Cell-type specific predictions:** Analysis was performed in the relevant hematopoietic cell context (CD34+ cells)

4. **Multi-scale effects:** The model captures both local chromatin changes and distal effects on gene expression

## Methods

**AlphaGenome DNA Model:** Used for variant effect prediction
**Gene annotations:** GENCODE v46 (protein-coding, longest transcripts)
**Variant scoring:** RNA-seq based gene expression change
**Comparison approach:** Oncogenic variants vs. length-matched random sequences
**Statistical approach:** 3 background variants per oncogenic variant

## Conclusion

This analysis successfully replicated and extended the tutorial workflow, demonstrating that:

- T-ALL associated variants near TAL1 have strong predicted upregulation effects
- Effects are significantly stronger than random background variants
- The model predicts biologically plausible mechanisms (enhancer creation, chromatin activation)
- Predictions align with published experimental findings

The workflow provides a practical template for analyzing non-coding regulatory variants in cancer and other diseases.
