#!/usr/bin/env python3
"""
Quick analysis of batch variant scoring results
"""

import pandas as pd
from pathlib import Path


def main():
    results_dir = Path('/home/kyuwon/projects/alphagenome/results/batch_variant_scoring')

    # Load high impact variants
    df = pd.read_csv(results_dir / 'high_impact_variants.csv')

    print("=" * 80)
    print("BATCH VARIANT SCORING - QUICK ANALYSIS")
    print("=" * 80)

    print(f"\n1. HIGH IMPACT VARIANTS SUMMARY")
    print(f"   Total high impact scores: {len(df)}")
    print(f"   Threshold: |raw_score| > 0.01")

    print(f"\n2. DISTRIBUTION BY VARIANT")
    variant_counts = df['variant_id'].value_counts()
    for variant, count in variant_counts.items():
        print(f"   {variant}: {count} high impact scores")

    print(f"\n3. DISTRIBUTION BY OUTPUT TYPE (assay)")
    output_counts = df['output_type'].value_counts()
    for output_type, count in output_counts.items():
        print(f"   {output_type}: {count} scores")

    print(f"\n4. TOP 10 POSITIVE EFFECTS (all variants)")
    top_positive = df.nlargest(10, 'raw_score')[
        ['variant_id', 'biosample_name', 'output_type', 'raw_score', 'quantile_score']
    ]
    for idx, row in top_positive.iterrows():
        print(f"   {row['variant_id']:20s} | {row['biosample_name']:30s} | "
              f"{row['output_type']:8s} | score: {row['raw_score']:7.4f} (q: {row['quantile_score']:6.3f})")

    print(f"\n5. TOP 10 NEGATIVE EFFECTS (all variants)")
    top_negative = df.nsmallest(10, 'raw_score')[
        ['variant_id', 'biosample_name', 'output_type', 'raw_score', 'quantile_score']
    ]
    for idx, row in top_negative.iterrows():
        print(f"   {row['variant_id']:20s} | {row['biosample_name']:30s} | "
              f"{row['output_type']:8s} | score: {row['raw_score']:7.4f} (q: {row['quantile_score']:6.3f})")

    print(f"\n6. MOST VARIABLE BIOSAMPLE TYPES")
    biosample_stats = df.groupby('biosample_type')['raw_score'].agg(['count', 'mean', 'std']).sort_values('std', ascending=False)
    print(biosample_stats.head(5))

    print(f"\n7. EFFECT SIZE DISTRIBUTION")
    print(f"   Positive scores (>0.01): {len(df[df['raw_score'] > 0.01])}")
    print(f"   Negative scores (<-0.01): {len(df[df['raw_score'] < -0.01])}")
    print(f"   Very high impact (|score| > 0.05): {len(df[abs(df['raw_score']) > 0.05])}")
    print(f"   Extreme impact (|score| > 0.1): {len(df[abs(df['raw_score']) > 0.1])}")

    print("\n" + "=" * 80)


if __name__ == '__main__':
    main()
