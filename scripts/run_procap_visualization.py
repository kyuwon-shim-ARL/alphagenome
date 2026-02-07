#!/usr/bin/env python3
"""
PROCAP Visualization Script
Generates RNA Polymerase II activity (PROCAP) predictions for 6 cell lines.

PROCAP (Precision Run-On Sequencing Cap) measures RNA Polymerase II activity
at transcription start sites with base-pair resolution.

Output:
    results/procap_visualization/
    ├── procap.png         # PROCAP visualization for all 6 cell lines
    └── results.json       # Metadata and results summary
"""

import json
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import gene_annotation, genome, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

load_dotenv()

# Configuration
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/procap_visualization")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

API_DELAY = 2.0

# PROCAP supports 6 cell lines (as per specification)
# Omitting ontology_terms will retrieve all available cell lines
PROCAP_CELL_LINES = {
    'A673': 'EFO:0002106',
    'Caco-2': 'EFO:0001099',
    'K562': 'EFO:0002067',
    'Calu3': 'EFO:0002819',
    'MCF 10A': 'EFO:0001200',
    'HUVEC': 'CL:0002618',
}

# Gene annotation URL
HG38_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)


def get_api_key():
    """Retrieve API key from environment"""
    api_key = os.getenv('ALPHAGENOME_API_KEY')
    if not api_key:
        raise ValueError(
            "ALPHAGENOME_API_KEY environment variable not set. "
            "Please set it before running this script."
        )
    return api_key


def save_plot(fig, name, results):
    """Save plot and record in results"""
    filepath = RESULTS_DIR / f"{name}.png"
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    results[name] = {
        'file': str(filepath),
        'status': 'success'
    }
    print(f"  Saved: {filepath}")
    return filepath


def main():
    results = {
        'plots': {},
        'metadata': {},
    }

    try:
        print("=" * 70)
        print("AlphaGenome PROCAP Visualization")
        print("=" * 70)

        # Initialize
        api_key = get_api_key()
        dna_model = dna_client.create(api_key)

        # Load gene annotations
        print("\nLoading gene annotations...")
        gtf = pd.read_feather(HG38_GTF_FEATHER)
        gtf_transcript = gene_annotation.filter_transcript_support_level(
            gene_annotation.filter_protein_coding(gtf), ['1']
        )
        transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)

        # Define genomic interval (chr22 around APOL4 gene - same as visualization_tour)
        interval = genome.Interval('chr22', 36_150_498, 36_252_898).resize(
            dna_client.SEQUENCE_LENGTH_1MB
        )
        print(f"\nWorking with interval: {interval}")

        # Get transcripts for annotation
        transcripts = transcript_extractor.extract(interval)
        apol_transcripts = [t for t in transcripts
                           if t.info.get('gene_name', '').startswith('APOL')]

        results['metadata']['interval'] = str(interval)
        results['metadata']['cell_lines'] = PROCAP_CELL_LINES

        # ==================== PROCAP Visualization ====================
        print("\nGenerating PROCAP visualization...")
        print("  Requesting RNA Polymerase II activity predictions...")

        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.PROCAP},
            ontology_terms=list(PROCAP_CELL_LINES.values()),
        )
        time.sleep(API_DELAY)

        if output.procap is not None and output.procap.values.shape[-1] > 0:
            # Record cell lines that were retrieved
            retrieved_cell_lines = output.procap.metadata['biosample_name'].unique().tolist()
            results['metadata']['retrieved_cell_lines'] = retrieved_cell_lines
            print(f"  Retrieved data for {len(retrieved_cell_lines)} cell lines:")
            for cell_line in retrieved_cell_lines:
                print(f"    - {cell_line}")

            # Record track count
            track_count = output.procap.values.shape[-1]
            results['metadata']['track_count'] = track_count
            print(f"  Total tracks: {track_count}")

            # Create visualization
            fig = plot_components.plot(
                [
                    plot_components.TranscriptAnnotation(apol_transcripts),
                    plot_components.Tracks(
                        tdata=output.procap,
                        ylabel_template='PROCAP: {biosample_name} ({strand})',
                    ),
                ],
                interval=interval,
                title='Predicted RNA Polymerase II Activity (PROCAP)',
            )
            save_plot(fig, 'procap', results['plots'])

        else:
            print("  Warning: No PROCAP data available")
            results['plots']['procap'] = {'status': 'no_data'}

        # Save results summary
        results['summary'] = {
            'total_plots': len([p for p in results['plots'].values() if p.get('status') == 'success']),
            'expected_cell_lines': list(PROCAP_CELL_LINES.keys()),
            'retrieved_cell_lines': results['metadata'].get('retrieved_cell_lines', []),
        }

        results_file = RESULTS_DIR / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        # Print summary
        print("\n" + "=" * 70)
        print("PROCAP VISUALIZATION COMPLETE")
        print("=" * 70)
        print(f"\nResults saved to: {RESULTS_DIR}")
        print(f"\nPlots generated:")
        for name, info in results['plots'].items():
            status = info.get('status', 'unknown')
            if status == 'success':
                print(f"  [OK] {name}.png")
            else:
                print(f"  [--] {name}: {status}")

        print(f"\nSummary: {results_file}")

        return 0

    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

        # Save partial results
        results_file = RESULTS_DIR / 'results.json'
        results['error'] = str(e)
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        return 1


if __name__ == '__main__':
    sys.exit(main())
