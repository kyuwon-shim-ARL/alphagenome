#!/usr/bin/env python3
"""
Batch Variant Scoring Script
Extracted from tutorials/batch_variant_scoring.ipynb
Scores multiple variants using AlphaGenome and saves results to CSV/JSON.
"""

import os
import json
from io import StringIO
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers


def main():
    # Setup
    print("Setting up AlphaGenome DNA model...")
    api_key = os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        raise ValueError("ALPHAGENOME_API_KEY environment variable not set")

    dna_model = dna_client.create(api_key)

    # Define variants (limiting to 5 for API management)
    vcf_file = """variant_id\tCHROM\tPOS\tREF\tALT
chr3_58394738_A_T_b38\tchr3\t58394738\tA\tT
chr8_28520_G_C_b38\tchr8\t28520\tG\tC
chr16_636337_G_A_b38\tchr16\t636337\tG\tA
chr16_1135446_G_T_b38\tchr16\t1135446\tG\tT
chr1_100000_C_G_b38\tchr1\t100000\tC\tG
"""

    print("\nLoading variants from VCF data...")
    vcf = pd.read_csv(StringIO(vcf_file), sep='\t')
    print(f"Loaded {len(vcf)} variants:")
    print(vcf)

    # Validate required columns
    required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
    for column in required_columns:
        if column not in vcf.columns:
            raise ValueError(f'VCF file is missing required column: {column}.')

    # Configure parameters
    organism = dna_client.Organism.HOMO_SAPIENS
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_1MB']

    # Select scorers for batch processing.
    # RECOMMENDED_VARIANT_SCORERS contains 19 keys:
    #   11 base: ATAC, CAGE, CHIP_HISTONE, CHIP_TF, CONTACT_MAPS, DNASE,
    #            PROCAP, RNA_SEQ, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS
    #   1 scorer-only (no corresponding OutputType): POLYADENYLATION
    #   7 _ACTIVE variants: ATAC_ACTIVE, CAGE_ACTIVE, CHIP_HISTONE_ACTIVE,
    #                       CHIP_TF_ACTIVE, DNASE_ACTIVE, PROCAP_ACTIVE, RNA_SEQ_ACTIVE
    #
    # Enabled (6 OutputType scorers): RNA_SEQ, CAGE, ATAC, DNASE, CHIP_HISTONE, SPLICE_SITES
    # Disabled: CHIP_TF, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, CONTACT_MAPS, PROCAP
    #   - These require specific ontology terms or produce large outputs
    #   - POLYADENYLATION disabled (scorer-only, no OutputType)
    #   - _ACTIVE variants auto-selected via case-insensitive key matching
    print("\nConfiguring variant scorers...")
    scorer_selections = {
        'rna_seq': True,
        'cage': True,
        'atac': True,
        'dnase': True,
        'chip_histone': True,
        'polyadenylation': False,
        'splice_sites': True,
        'splice_site_usage': False,
        'splice_junctions': False,
    }

    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_scorers = [
        all_scorers[key]
        for key in all_scorers
        if scorer_selections.get(key.lower(), False)
    ]

    # Remove unsupported scorers
    unsupported_scorers = [
        scorer
        for scorer in selected_scorers
        if (
            organism.value
            not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
        )
    ]

    if len(unsupported_scorers) > 0:
        print(f'Excluding {unsupported_scorers} scorers as they are not supported for {organism}.')
        for unsupported_scorer in unsupported_scorers:
            selected_scorers.remove(unsupported_scorer)

    print(f"Using {len(selected_scorers)} scorers: {[s.base_variant_scorer.name for s in selected_scorers]}")

    # Score variants
    print("\nScoring variants (this may take a few minutes)...")
    results = []

    for i, vcf_row in tqdm(vcf.iterrows(), total=len(vcf), desc="Scoring variants"):
        variant = genome.Variant(
            chromosome=str(vcf_row.CHROM),
            position=int(vcf_row.POS),
            reference_bases=vcf_row.REF,
            alternate_bases=vcf_row.ALT,
            name=vcf_row.variant_id,
        )
        interval = variant.reference_interval.resize(sequence_length)

        variant_scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=selected_scorers,
            organism=organism,
        )
        results.append(variant_scores)

    # Convert results to DataFrame
    print("\nProcessing results...")
    df_scores = variant_scorers.tidy_scores(results)
    print(f"Generated {len(df_scores)} score rows")

    # Display summary statistics
    print("\nSummary Statistics:")
    print(f"Total variants scored: {vcf['variant_id'].nunique()}")
    print(f"Total score rows: {len(df_scores)}")
    print(f"Score columns: {df_scores.columns.tolist()}")

    # Check for cell_type or ontology_curie
    if 'cell_type' in df_scores.columns:
        print(f"Unique cell types: {df_scores['cell_type'].nunique()}")
    elif 'ontology_curie' in df_scores.columns:
        print(f"Unique ontology curies: {df_scores['ontology_curie'].nunique()}")

    # Show raw_score distribution
    print("\nRaw Score Statistics:")
    print(df_scores['raw_score'].describe())

    # Save results
    results_dir = Path('/home/kyuwon/projects/alphagenome/results/batch_variant_scoring')
    results_dir.mkdir(parents=True, exist_ok=True)

    csv_path = results_dir / 'variant_scores.csv'
    print(f"\nSaving results to CSV: {csv_path}")
    df_scores.to_csv(csv_path, index=False)

    # Create summary JSON with aggregated scores per variant
    print("\nCreating summary JSON...")
    summary_data = []

    # Determine which column to use for cell identification
    cell_col = 'cell_type' if 'cell_type' in df_scores.columns else 'ontology_curie'

    for variant_id in vcf['variant_id']:
        variant_df = df_scores[df_scores['variant_id'] == variant_id]

        if len(variant_df) == 0:
            continue

        # Build top scores columns dynamically
        score_cols = [col for col in [cell_col, 'variant_scorer', 'raw_score', 'quantile_score']
                      if col in variant_df.columns]

        summary_data.append({
            'variant_id': variant_id,
            'chromosome': variant_df['scored_interval'].iloc[0].split(':')[0] if len(variant_df) > 0 else None,
            'total_scores': len(variant_df),
            'mean_raw_score': float(variant_df['raw_score'].mean()),
            'median_raw_score': float(variant_df['raw_score'].median()),
            'std_raw_score': float(variant_df['raw_score'].std()),
            'mean_quantile_score': float(variant_df['quantile_score'].mean()) if 'quantile_score' in variant_df.columns else None,
            f'unique_{cell_col}s': int(variant_df[cell_col].nunique()) if cell_col in variant_df.columns else 0,
            'unique_scorers': int(variant_df['variant_scorer'].nunique()) if 'variant_scorer' in variant_df.columns else 0,
            'top_positive_scores': variant_df.nlargest(3, 'raw_score')[score_cols].to_dict('records'),
            'top_negative_scores': variant_df.nsmallest(3, 'raw_score')[score_cols].to_dict('records'),
        })

    json_path = results_dir / 'variant_scores_summary.json'
    print(f"Saving summary to JSON: {json_path}")
    print(f"Summary data entries: {len(summary_data)}")
    with open(json_path, 'w') as f:
        json.dump(summary_data, f, indent=2, default=str)

    # Create filtered view: high impact variants (|raw_score| > 0.01)
    high_impact_df = df_scores[abs(df_scores['raw_score']) > 0.01]
    high_impact_path = results_dir / 'high_impact_variants.csv'
    print(f"\nSaving high impact variants (|raw_score| > 0.01): {high_impact_path}")
    print(f"High impact rows: {len(high_impact_df)}")
    high_impact_df.to_csv(high_impact_path, index=False)

    # Display sample of results
    print("\nSample results (first 10 rows):")
    print(df_scores.head(10).to_string())

    print("\n" + "="*80)
    print("BATCH VARIANT SCORING COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {results_dir}")
    print(f"  - Full scores: {csv_path}")
    print(f"  - Summary: {json_path}")
    print(f"  - High impact: {high_impact_path}")


if __name__ == '__main__':
    main()
