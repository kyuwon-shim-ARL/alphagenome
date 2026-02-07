#!/usr/bin/env python3
"""
AlphaGenome ISM Analysis - 256bp Extension
Extends In Silico Mutagenesis analysis from 64bp to 256bp
"""

import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.interpretation import ism
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components

# Load environment variables from .env file
load_dotenv()

# Configuration
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/ism_256bp")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# API delay to avoid rate limiting
API_DELAY = 2.0


def get_api_key():
    """Retrieve API key from environment"""
    api_key = os.getenv('ALPHAGENOME_API_KEY')
    if not api_key:
        raise ValueError(
            "ALPHAGENOME_API_KEY environment variable not set. "
            "Please set it before running this script."
        )
    return api_key


def main():
    results = {
        "analysis": "ism_256bp",
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "ism_width": 256,
            "expected_variants": 768,  # 256bp x 3 alternate alleles
            "api_delay": API_DELAY,
        }
    }

    try:
        print("="*60)
        print("AlphaGenome ISM Analysis - 256bp Extension")
        print("="*60)

        print("\nInitializing AlphaGenome client...")
        api_key = get_api_key()
        dna_model = dna_client.create(api_key)
        print("✓ Client initialized")

        # Define genomic interval (same region as quick_start example)
        print("\nDefining genomic interval...")
        sequence_interval = genome.Interval('chr20', 3_753_000, 3_753_400)
        sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
        print(f"✓ Sequence interval: {sequence_interval}")

        # ISM interval - 256bp (extended from 64bp)
        ism_interval = sequence_interval.resize(256)
        print(f"✓ ISM interval: {ism_interval}")
        print(f"  Width: {ism_interval.width} bp")
        print(f"  Expected variants: {ism_interval.width * 3} (3 alternate alleles per position)")

        results["interval"] = {
            "sequence_interval": str(sequence_interval),
            "ism_interval": str(ism_interval),
            "ism_width": ism_interval.width,
        }

        # Create CenterMaskScorer for DNASE
        print("\nConfiguring variant scorer...")
        dnase_variant_scorer = variant_scorers.CenterMaskScorer(
            requested_output=dna_client.OutputType.DNASE,
            width=501,
            aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
        )
        print("✓ CenterMaskScorer configured for DNASE")

        # Score ISM variants (batch processed)
        print(f"\nScoring ISM variants...")
        print(f"  This may take several minutes due to {ism_interval.width * 3} variants...")
        print(f"  API delay: {API_DELAY}s between requests")

        start_time = time.time()
        ism_scores = dna_model.score_ism_variants(
            interval=sequence_interval,
            ism_interval=ism_interval,
            variant_scorers=[dnase_variant_scorer],
        )
        elapsed_time = time.time() - start_time

        print(f"✓ ISM scoring completed in {elapsed_time:.1f}s")
        print(f"  Variants scored: {len(ism_scores)}")

        results["scoring"] = {
            "variants_scored": len(ism_scores),
            "elapsed_time_seconds": elapsed_time,
        }

        # Convert ISM scores to DataFrame format (manual conversion - tidy_scores incompatible with ISM)
        print("\nConverting ISM scores to DataFrame...")
        # ISM scores structure: [(AnnData,), (AnnData,), ...] - one tuple per variant
        # Each AnnData has .X (scores), .uns['variant'] (Variant object)
        ism_rows = []
        for score_tuple in ism_scores:
            adata = score_tuple[0]  # Get first (and only) scorer's AnnData
            variant = adata.uns['variant']
            mean_score = float(adata.X.mean())
            ism_rows.append({
                'chromosome': variant.chromosome,
                'position': variant.position,
                'reference': variant.reference_bases,
                'alternate': variant.alternate_bases,
                'raw_score': mean_score,
            })
        ism_df = pd.DataFrame(ism_rows)
        print(f"✓ ISM scores DataFrame shape: {ism_df.shape}")
        print(f"  Columns: {list(ism_df.columns)}")

        results["dataframe"] = {
            "shape": list(ism_df.shape),
            "columns": list(ism_df.columns),
        }

        # Save ISM scores to CSV
        csv_path = RESULTS_DIR / "ism_scores.csv"
        ism_df.to_csv(csv_path, index=False)
        print(f"✓ Saved ISM scores: {csv_path}")

        # Create ISM matrix for visualization
        print("\nGenerating ISM matrix...")
        # Extract mean score across all DNASE tracks for each variant
        # ism_scores structure: [(AnnData, ), (AnnData, ), ...] - one tuple per variant
        def extract_mean_score(score_tuple):
            """Extract mean score from AnnData object"""
            adata = score_tuple[0]  # Get first (and only) scorer's AnnData
            return float(adata.X.mean())  # Mean across all tracks

        ism_matrix_data = ism.ism_matrix(
            variant_scores=[extract_mean_score(score_tuple) for score_tuple in ism_scores],
            variants=[score_tuple[0].uns['variant'] for score_tuple in ism_scores],
        )
        print(f"✓ ISM matrix shape: {ism_matrix_data.shape}")

        # Save ISM heatmap using SeqLogo
        print("\nGenerating ISM heatmap...")
        plot_components.plot(
            [
                plot_components.SeqLogo(
                    scores=ism_matrix_data,
                    scores_interval=ism_interval,
                    ylabel='ISM DNASE (mean)',
                )
            ],
            interval=ism_interval,
            fig_width=35,
        )
        heatmap_path = RESULTS_DIR / "ism_heatmap.png"
        plt.savefig(heatmap_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved ISM heatmap: {heatmap_path}")

        results["outputs"] = {
            "csv": str(csv_path),
            "heatmap": str(heatmap_path),
        }

        # Compute summary statistics
        print("\nComputing summary statistics...")
        # CenterMaskScorer doesn't have genes, so we just compute stats on raw_score
        if 'raw_score' in ism_df.columns:
            score_stats = {
                "mean": float(ism_df['raw_score'].mean()),
                "std": float(ism_df['raw_score'].std()),
                "min": float(ism_df['raw_score'].min()),
                "max": float(ism_df['raw_score'].max()),
                "median": float(ism_df['raw_score'].median()),
            }
            print(f"  Mean raw_score: {score_stats['mean']:.6f}")
            print(f"  Std raw_score: {score_stats['std']:.6f}")
            print(f"  Min raw_score: {score_stats['min']:.6f}")
            print(f"  Max raw_score: {score_stats['max']:.6f}")
            print(f"  Median raw_score: {score_stats['median']:.6f}")
        else:
            print(f"  Warning: 'raw_score' column not found in DataFrame")
            print(f"  Available columns: {list(ism_df.columns)}")
            score_stats = {}

        results["statistics"] = score_stats

        # Save results summary
        results["status"] = "success"
        results_file = RESULTS_DIR / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)

        print(f"\n{'='*60}")
        print(f"✓ ISM 256bp analysis completed successfully!")
        print(f"✓ Results saved to: {RESULTS_DIR}")
        print(f"✓ Summary: {results_file}")
        print(f"{'='*60}")

        # Display key result
        print(f"\nKey Result: ISM Variants Scored = {len(ism_scores)}")

        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

        # Save error information
        results["status"] = "error"
        results["error"] = str(e)
        results_file = RESULTS_DIR / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)

        return 1


if __name__ == '__main__':
    sys.exit(main())
