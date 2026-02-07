# Batch Variant Scoring Results

## Overview

This directory contains the results from running the AlphaGenome batch variant scoring tutorial. The analysis scored 5 genetic variants across multiple cell types and tissue contexts using AlphaGenome's DNA foundation model.

## Execution Summary

- **Date**: 2026-02-05
- **Tutorial**: `tutorials/batch_variant_scoring.ipynb`
- **Script**: `scripts/run_batch_variant_scoring.py`
- **Variants Analyzed**: 5
- **Total Score Rows**: 121,550
- **Unique Cell Types/Contexts**: 687 (ontology curies)
- **Scorers Used**: 5 (ATAC, CAGE, DNASE, ChIP-histone, RNA-seq)

## Variants Analyzed

1. **chr3:58394738:A>T** (chr3_58394738_A_T_b38)
2. **chr8:28520:G>C** (chr8_28520_G_C_b38)
3. **chr16:636337:G>A** (chr16_636337_G_A_b38)
4. **chr16:1135446:G>T** (chr16_1135446_G_T_b38)
5. **chr1:100000:C>G** (chr1_100000_C_G_b38)

## Output Files

### 1. `variant_scores.csv` (30 MB, 121,550 rows)

Complete scoring results with columns:
- `variant_id`: Variant identifier
- `scored_interval`: Genomic interval used for scoring
- `output_type`: Type of assay (ATAC, CAGE, DNASE, etc.)
- `variant_scorer`: Scorer algorithm and parameters
- `ontology_curie`: Cell/tissue ontology identifier (CL, EFO, UBERON)
- `biosample_name`: Human-readable cell/tissue name
- `biosample_type`: Type (primary_cell, cell_line, tissue, etc.)
- `raw_score`: Raw variant effect score
- `quantile_score`: Normalized quantile score (-1 to 1)

### 2. `high_impact_variants.csv` (1.1 MB, 4,509 rows)

Filtered results showing only high-impact scores where |raw_score| > 0.01.

### 3. `variant_scores_summary.json` (2 KB)

Aggregated statistics per variant including:
- Total number of scores
- Mean, median, and standard deviation of raw scores
- Top 3 positive and negative scores
- Unique cell types and scorers used

## Key Findings

### Variant Impact Summary

| Variant ID | Total Scores | Mean Raw Score | Median Raw Score | Std Dev |
|------------|--------------|----------------|------------------|---------|
| chr3:58394738:A>T | 12,430 | 0.000440 | 0.000117 | 0.004936 |
| chr8:28520:G>C | 10,450 | 0.000029 | 0.000010 | 0.001697 |
| chr16:636337:G>A | 40,150 | 0.003363 | -0.000000 | 0.018540 |
| chr16:1135446:G>T | 41,734 | -0.000841 | -0.000020 | 0.010343 |
| chr1:100000:C>G | 16,786 | 0.000676 | 0.000007 | 0.004511 |

### Overall Statistics

- **Mean raw score**: 0.000963
- **Median raw score**: 0.000000
- **Standard deviation**: 0.012607
- **Min score**: -0.653481
- **Max score**: 0.313484
- **High impact scores** (|score| > 0.01): 4,509 (3.7% of total)

### Most Impactful Variant: chr16:636337:G>A

This variant showed:
- Highest mean raw score (0.003363)
- Highest standard deviation (0.018540)
- Most variable effects across cell types

### Top Positive Effects (chr3:58394738:A>T)

1. **UBERON:0000955** (DNase): +0.0560, quantile: 0.799
2. **EFO:0007950 (GM23338 cell line)** (ATAC): +0.0472, quantile: 0.768
3. **EFO:0007598** (DNase): +0.0457, quantile: 0.840

### Top Negative Effects (chr3:58394738:A>T)

1. **CL:0001059** (DNase): -0.1229, quantile: -0.927
2. **EFO:0003037** (DNase): -0.0912, quantile: -0.920
3. **EFO:0002067 (K562 cell line)** (DNase): -0.0692, quantile: -0.892

## Methodology

### Sequence Length
- 1 MB context window around each variant

### Variant Scorers Used
1. **CENTER_MASK (ATAC-seq)**: Chromatin accessibility changes
2. **CENTER_MASK (CAGE)**: Transcription start site effects
3. **CENTER_MASK (DNASE)**: DNase hypersensitivity changes
4. **CENTER_MASK (ChIP-histone)**: Histone modification effects
5. **GENE_MASK_LFC (RNA-seq)**: Gene expression log-fold changes

### Scoring Parameters
- **Aggregation**: DIFF_LOG2_SUM (difference between reference and alternate)
- **Width**: 501 bp center mask
- **Organism**: Homo sapiens (hg38)

## Usage Examples

### Loading Results in Python

```python
import pandas as pd

# Load full results
df = pd.read_csv('variant_scores.csv')

# Filter to specific cell type
t_cells = df[df['ontology_curie'] == 'CL:0000084']

# Get high impact scores
high_impact = df[abs(df['raw_score']) > 0.01]

# Analyze by variant
variant_summary = df.groupby('variant_id')['raw_score'].describe()
```

### Analyzing Specific Variants

```python
import json

with open('variant_scores_summary.json') as f:
    summary = json.load(f)

for variant in summary:
    print(f"{variant['variant_id']}: mean={variant['mean_raw_score']:.6f}")
```

## Notes

- Variants were limited to 5 for API management
- Total API calls: 5 variants × ~24,000 scores each = ~7 seconds
- All scores use hg38 reference genome
- Ontology identifiers: CL (Cell), EFO (Experimental Factor), UBERON (Anatomy)

## Related Files

- **Tutorial notebook**: `/home/kyuwon/projects/alphagenome/tutorials/batch_variant_scoring.ipynb`
- **Execution script**: `/home/kyuwon/projects/alphagenome/scripts/run_batch_variant_scoring.py`
- **Summary generator**: `/home/kyuwon/projects/alphagenome/scripts/generate_batch_summary.py`
