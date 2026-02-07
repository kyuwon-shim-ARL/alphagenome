# AlphaGenome Visualization Modality Tour Results

## Overview

This directory contains visualization outputs from the AlphaGenome Visualization Modality Tour tutorial, demonstrating the model's ability to predict various genomic data modalities.

## Genomic Region

- **Chromosome**: chr22
- **Interval**: 35,677,410 - 36,725,986 (1,048,576 bp / 1 Mb)
- **Primary Genes**: APOL1, APOL2, APOL3, APOL4, APOL5, APOL6 (Apolipoprotein family)

## Generated Visualizations

### 1. RNA_SEQ - Gene Expression (01_rna_seq.png)
- **Size**: 118 KB
- **Dimensions**: 2722 x 817 px
- **Description**: RNA-seq predictions showing gene expression levels in colon tissue
- **Tissues**: Colon - Sigmoid, Colon - Transverse
- **Key Features**: Displays strand-specific expression patterns across the APOL gene family

### 2. CAGE - Transcription Start Sites (02_cage.png)
- **Size**: 59 KB
- **Dimensions**: 2590 x 470 px
- **Description**: CAGE predictions identifying transcription start sites
- **Tissues**: Colon - Sigmoid, Colon - Transverse
- **Key Features**: Sharp peaks at TSS locations, particularly around APOL1

### 3. DNASE - Chromatin Accessibility (03_dnase.png)
- **Size**: 60 KB
- **Dimensions**: 2666 x 355 px
- **Description**: DNase-seq predictions showing open chromatin regions
- **Tissues**: Colon - Transverse, Colon - Sigmoid
- **Key Features**: Indicates regulatory regions and accessible DNA

### 4. ATAC - Open Chromatin (04_atac.png)
- **Size**: 58 KB
- **Dimensions**: Similar to DNASE
- **Description**: ATAC-seq predictions showing chromatin accessibility
- **Tissues**: Colon - Transverse, Colon - Sigmoid
- **Key Features**: Complementary to DNASE, shows accessible chromatin regions

### 5. CHIP_HISTONE - Histone Modifications (05_chip_histone.png)
- **Size**: 553 KB (largest file)
- **Dimensions**: Complex multi-track visualization
- **Description**: ChIP-seq predictions for histone modification marks
- **Tissues**: Multiple colon tissues
- **Histone Marks**:
  - H3K27AC (red) - Active enhancers
  - H3K36ME3 (orange) - Gene bodies
  - H3K4ME1 (blue) - Enhancers
  - H3K4ME3 (purple) - Active promoters
  - H3K9AC (green) - Active transcription
  - H3K27ME3 (pink) - Repressed regions
- **Key Features**: Color-coded by histone mark type, shows epigenetic landscape

### 6. SPLICE_SITES and SPLICE_JUNCTIONS (06_splice.png)
- **Size**: 80 KB
- **Dimensions**: 2722 x 817 px
- **Description**: Splicing predictions showing splice sites and junctions
- **Focus**: APOL4 gene region
- **Tissues**: Colon - Transverse, Colon - Sigmoid
- **Key Features**:
  - Splice sites (acceptor/donor) predictions
  - Splice junction arcs (Sashimi plot style)
  - Strand-specific visualization (negative strand for APOL4)

### 7. CONTACT_MAPS - 3D Chromatin Interaction (07_contact_maps.png)
- **Size**: 508 KB (second largest)
- **Dimensions**: Complex heatmap visualization
- **Description**: Hi-C contact map predictions showing 3D chromatin structure
- **Cell Line**: HCT116 (colon carcinoma)
- **Resolution**: 2048 bp
- **Key Features**: Topologically-associated domains (TADs) visible as blocks

## Technical Details

### API Configuration
- Model: AlphaGenome DNA Client
- Organism: Homo sapiens (human)
- Reference: hg38
- Gene Annotations: GENCODE v46
- Protein-coding genes with transcript support level 1

### Data Processing
- API delay between requests: 2 seconds
- Image format: PNG (150 DPI)
- Backend: matplotlib (Agg)

### Execution
- Script: `/home/kyuwon/projects/alphagenome/scripts/run_visualization_tour.py`
- Duration: ~1 minute
- Status: All 7 visualizations completed successfully

## Results Summary

```json
{
  "total": 7,
  "successful": 7,
  "genomic_interval": "chr22:35677410-36725986",
  "output_types": [
    "RNA_SEQ",
    "CAGE",
    "DNASE",
    "ATAC",
    "CHIP_HISTONE",
    "SPLICE_SITES/SPLICE_JUNCTIONS",
    "CONTACT_MAPS"
  ]
}
```

## Key Insights

1. **Gene Expression**: Strong expression signals for APOL1 in colon tissue (positive strand)
2. **Chromatin Accessibility**: Multiple accessible regions indicating regulatory elements
3. **Histone Modifications**: Promoter marks (H3K4me3) align with transcription start sites
4. **Splicing**: Complex splicing patterns in APOL4 gene with multiple exons
5. **3D Structure**: Visible TAD boundaries and chromatin interaction blocks

## Files Generated

- `01_rna_seq.png` - RNA expression visualization
- `02_cage.png` - CAGE transcription start sites
- `03_dnase.png` - DNase chromatin accessibility
- `04_atac.png` - ATAC chromatin accessibility
- `05_chip_histone.png` - Histone modification landscape
- `06_splice.png` - Splicing events
- `07_contact_maps.png` - 3D chromatin contacts
- `results.json` - Structured results summary
- `README.md` - This documentation file

## References

- Tutorial Source: `/home/kyuwon/projects/alphagenome/tutorials/visualization_modality_tour.ipynb`
- AlphaGenome Documentation: https://www.alphagenomedocs.com/
- GENCODE Annotations: https://www.gencodegenes.org/
