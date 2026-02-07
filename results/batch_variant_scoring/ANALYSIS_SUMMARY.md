# Batch Variant Scoring Analysis Summary

## Execution Details

**Date**: 2026-02-05
**Runtime**: ~7 seconds for 5 variants
**Total Scores Generated**: 121,550
**High Impact Scores**: 4,509 (3.7%)

## Key Findings

### 1. Most Impactful Variant: chr16:636337:G>A

This variant demonstrated the strongest effects:
- **2,424 high impact scores** (53.7% of all high impact)
- Strongest effect on immune cells, particularly:
  - Neutrophils (CAGE score: +0.3135, 99.5th percentile)
  - Eosinophils (CAGE score: +0.3001, 99.2nd percentile)
  - Basophils and Langerhans cells

### 2. Most Detrimental Variant: chr16:1135446:G>T

This variant showed strong negative regulatory effects:
- **1,133 high impact scores**
- Primarily affects myeloid cells:
  - CD14+ monocytes (ChIP-histone: -0.6535, -99.7th percentile)
  - Inflammatory macrophages (DNase: -0.5527)
  - Common myeloid progenitors (DNase: -0.5501)

### 3. Distribution of Effects

**By Assay Type**:
- ChIP-Histone: 1,718 high impact scores (38%)
- CAGE: 962 scores (21%)
- DNase: 807 scores (18%)
- RNA-seq: 684 scores (15%)
- ATAC: 338 scores (7%)

**By Biosample Type**:
- Primary cells: 1,300 scores (most variable, std=0.074)
- Tissues: 1,681 scores (std=0.055)
- Cell lines: 1,246 scores (std=0.052)
- In vitro differentiated: 282 scores (std=0.039)

### 4. Effect Magnitude Classification

| Category | Threshold | Count | Percentage |
|----------|-----------|-------|------------|
| High impact | \|score\| > 0.01 | 4,509 | 3.7% |
| Very high impact | \|score\| > 0.05 | 1,423 | 1.2% |
| Extreme impact | \|score\| > 0.1 | 546 | 0.4% |

**Directional Bias**:
- Positive effects: 3,079 (68.3%)
- Negative effects: 1,430 (31.7%)

## Top 10 Strongest Positive Effects

| Variant | Cell/Tissue | Assay | Raw Score | Quantile |
|---------|-------------|-------|-----------|----------|
| chr16:636337:G>A | Neutrophil | CAGE | +0.3135 | 99.5% |
| chr16:636337:G>A | Eosinophil | CAGE | +0.3001 | 99.2% |
| chr16:636337:G>A | Langerhans cell | CAGE | +0.2743 | 99.4% |
| chr16:636337:G>A | Basophil | CAGE | +0.2713 | 98.9% |
| chr16:636337:G>A | Neutrophil | CAGE | +0.2646 | 99.5% |
| chr16:636337:G>A | Knee ligament | CAGE | +0.2529 | 99.7% |
| chr16:636337:G>A | Basophil | CAGE | +0.2527 | 98.8% |
| chr16:636337:G>A | Calcaneal tendon | CAGE | +0.2466 | 99.8% |
| chr16:636337:G>A | Classical monocyte | CAGE | +0.2366 | 99.2% |
| chr16:636337:G>A | CD14+ monocyte | DNase | +0.2287 | 98.2% |

**Interpretation**: The chr16:636337:G>A variant shows exceptionally strong positive regulatory effects on CAGE (transcription start sites) in immune cells, particularly myeloid lineage cells. This suggests enhanced transcriptional activity.

## Top 10 Strongest Negative Effects

| Variant | Cell/Tissue | Assay | Raw Score | Quantile |
|---------|-------------|-------|-----------|----------|
| chr16:1135446:G>T | CD14+ monocyte | ChIP-Histone | -0.6535 | -99.7% |
| chr16:1135446:G>T | Inflammatory macrophage | DNase | -0.5527 | -99.2% |
| chr16:1135446:G>T | Myeloid progenitor | DNase | -0.5501 | -99.1% |
| chr16:1135446:G>T | Suppressor macrophage | DNase | -0.5208 | -99.3% |
| chr16:1135446:G>T | NB4 (cell line) | DNase | -0.4718 | -99.2% |
| chr16:1135446:G>T | HL-60 (cell line) | DNase | -0.4511 | -99.4% |
| chr16:1135446:G>T | KBM-7 (cell line) | DNase | -0.4213 | -99.1% |
| chr16:1135446:G>T | Myeloid progenitor | ChIP-Histone | -0.3764 | -99.6% |
| chr16:1135446:G>T | CD14+ monocyte | ChIP-Histone | -0.3733 | -99.5% |
| chr16:1135446:G>T | CD14+ monocyte | DNase | -0.3486 | -99.1% |

**Interpretation**: The chr16:1135446:G>T variant disrupts regulatory elements in myeloid cells, particularly affecting chromatin accessibility (DNase) and histone modifications. This suggests reduced regulatory activity and potential gene silencing.

## Biological Insights

### 1. Cell Type Specificity

Both high-impact variants (chr16:636337:G>A and chr16:1135446:G>T) show strong myeloid cell-specific effects, suggesting these variants may:
- Affect myeloid differentiation pathways
- Influence immune cell function
- Impact inflammatory responses

### 2. Regulatory Mechanisms

The dominant assay types affected are:
- **CAGE (21%)**: Transcription start site changes
- **ChIP-Histone (38%)**: Histone modification alterations
- **DNase (18%)**: Chromatin accessibility changes

This indicates the variants primarily affect transcriptional regulation rather than RNA stability or splicing.

### 3. Tissue Context

Primary cells show the most variable responses (std=0.074), suggesting:
- Context-dependent regulatory effects
- Potential environmental modifiers
- Cell state-specific impacts

## Recommendations for Follow-up

1. **Functional validation**: Prioritize chr16:636337:G>A and chr16:1135446:G>T for experimental validation in myeloid cells

2. **GWAS overlap**: Check if these variants are near known immune disease risk loci

3. **Gene annotation**: Identify nearest genes and regulatory elements at chr16:636337 and chr16:1135446

4. **Cell-type experiments**: Focus on monocytes, neutrophils, and macrophages for functional assays

5. **Mechanism studies**: Investigate CAGE-seq and ChIP-seq changes in relevant cell types

## Technical Notes

- All scores normalized to quantiles based on genome-wide variant effects
- 1 MB sequence context used for prediction
- hg38 reference genome
- Multiple independent scorers provide orthogonal evidence

## Files Generated

1. `variant_scores.csv` - Complete results (121,550 scores)
2. `high_impact_variants.csv` - Filtered high-impact scores (4,509 scores)
3. `variant_scores_summary.json` - Aggregated statistics per variant
4. `README.md` - Methodology and overview
5. `ANALYSIS_SUMMARY.md` - This file

## Citation

Data generated using AlphaGenome DNA Foundation Model (Google DeepMind, 2025).
