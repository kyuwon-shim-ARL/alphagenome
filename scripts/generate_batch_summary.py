#!/usr/bin/env python3
"""
Generate summary JSON from batch variant scoring CSV results
"""

import json
from pathlib import Path
import pandas as pd


def main():
    results_dir = Path('/home/kyuwon/projects/alphagenome/results/batch_variant_scoring')
    csv_path = results_dir / 'variant_scores.csv'

    print(f"Loading results from: {csv_path}")
    df_scores = pd.read_csv(csv_path)

    print(f"Loaded {len(df_scores)} rows")
    print(f"Columns: {df_scores.columns.tolist()}")

    # Get unique variants
    variants = df_scores['variant_id'].unique()
    print(f"Unique variants: {len(variants)}")

    # Determine cell identifier column
    cell_col = 'cell_type' if 'cell_type' in df_scores.columns else 'ontology_curie'

    summary_data = []

    for variant_id in variants:
        variant_df = df_scores[df_scores['variant_id'] == variant_id]

        if len(variant_df) == 0:
            continue

        # Build score columns dynamically
        score_cols = [col for col in [cell_col, 'variant_scorer', 'raw_score', 'quantile_score']
                      if col in variant_df.columns]

        summary_entry = {
            'variant_id': variant_id,
            'chromosome': variant_df['scored_interval'].iloc[0].split(':')[0],
            'total_scores': int(len(variant_df)),
            'mean_raw_score': float(variant_df['raw_score'].mean()),
            'median_raw_score': float(variant_df['raw_score'].median()),
            'std_raw_score': float(variant_df['raw_score'].std()),
        }

        if 'quantile_score' in variant_df.columns:
            summary_entry['mean_quantile_score'] = float(variant_df['quantile_score'].mean())

        if cell_col in variant_df.columns:
            summary_entry[f'unique_{cell_col}s'] = int(variant_df[cell_col].nunique())

        if 'variant_scorer' in variant_df.columns:
            summary_entry['unique_scorers'] = int(variant_df['variant_scorer'].nunique())

        # Top positive and negative scores
        top_pos = variant_df.nlargest(3, 'raw_score')[score_cols]
        top_neg = variant_df.nsmallest(3, 'raw_score')[score_cols]

        summary_entry['top_positive_scores'] = top_pos.to_dict('records')
        summary_entry['top_negative_scores'] = top_neg.to_dict('records')

        summary_data.append(summary_entry)

    # Save summary
    json_path = results_dir / 'variant_scores_summary.json'
    print(f"\nSaving {len(summary_data)} variant summaries to: {json_path}")

    with open(json_path, 'w') as f:
        json.dump(summary_data, f, indent=2, default=str)

    print("\nSummary Statistics:")
    for entry in summary_data:
        print(f"\n{entry['variant_id']}:")
        print(f"  Total scores: {entry['total_scores']}")
        print(f"  Mean raw score: {entry['mean_raw_score']:.6f}")
        print(f"  Median raw score: {entry['median_raw_score']:.6f}")
        print(f"  Std raw score: {entry['std_raw_score']:.6f}")


if __name__ == '__main__':
    main()
