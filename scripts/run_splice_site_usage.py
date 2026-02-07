#!/usr/bin/env python3
"""
SPLICE_SITE_USAGE Standalone Visualization Script
Generates SPLICE_SITE_USAGE visualization for APOL4 gene region
"""

import json
import os
import sys
import time
import traceback
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import gene_annotation, genome, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

# Load environment variables from .env file
load_dotenv()

# Configuration
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/splice_site_usage")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# API delay to avoid rate limiting
API_DELAY = 2.0


def get_api_key():
    """Retrieve API key from environment"""
    api_key = os.getenv('ALPHAGENOME_API_KEY')
    if not api_key:
        raise ValueError(
            "ALPHAGENOME_API_KEY environment variable not set. "
            "Please set it in .env file before running this script."
        )
    return api_key


def main():
    results = {
        'interval': None,
        'gene': 'APOL4',
        'ontology_terms': [],
        'outputs': {},
        'status': 'running'
    }

    try:
        print("="*60)
        print("SPLICE_SITE_USAGE Visualization Script")
        print("="*60)

        # Initialize client
        print("\n[1/4] Initializing AlphaGenome client...")
        api_key = get_api_key()
        dna_model = dna_client.create(api_key)
        print("✓ Client initialized")

        # Load gene annotations
        print("\n[2/4] Loading gene annotations...")
        gtf = pd.read_feather(
            'https://storage.googleapis.com/alphagenome/reference/gencode/'
            'hg38/gencode.v46.annotation.gtf.gz.feather'
        )

        # Filter to protein-coding genes and highly supported transcripts
        gtf_transcript = gene_annotation.filter_transcript_support_level(
            gene_annotation.filter_protein_coding(gtf), ['1']
        )

        # Create transcript extractor
        transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
        print("✓ Gene annotations loaded")

        # Define genomic interval (chr22 around APOL4 gene)
        interval = genome.Interval('chr22', 36_150_498, 36_252_898).resize(
            dna_client.SEQUENCE_LENGTH_1MB
        )
        results['interval'] = str(interval)
        print(f"✓ Working with interval: {interval}")

        # Define ontology terms (Colon Transverse, Colon Sigmoid)
        ontology_terms = [
            'UBERON:0001157',  # Colon - Transverse
            'UBERON:0001159',  # Colon - Sigmoid
        ]
        results['ontology_terms'] = ontology_terms

        # Make predictions
        print("\n[3/4] Making predictions...")
        print(f"  Requesting: SPLICE_SITES, SPLICE_SITE_USAGE")
        print(f"  Ontology terms: {ontology_terms}")

        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={
                dna_client.OutputType.SPLICE_SITES,
                dna_client.OutputType.SPLICE_SITE_USAGE,
            },
            ontology_terms=ontology_terms,
        )
        time.sleep(API_DELAY)
        print("✓ Predictions completed")

        # Get APOL4 interval
        print("\n[4/4] Generating visualization...")
        apol4_interval = gene_annotation.get_gene_interval(gtf, gene_symbol='APOL4')
        apol4_interval.resize_inplace(apol4_interval.width + 1000)
        print(f"  Focus interval: {apol4_interval}")

        # Extract all transcripts in the region
        all_transcripts = transcript_extractor.extract(interval)

        # Build plot with TranscriptAnnotation, SPLICE_SITES, and SPLICE_SITE_USAGE
        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(all_transcripts),
                plot_components.Tracks(
                    tdata=output.splice_sites.filter_to_negative_strand(),
                    ylabel_template='SPLICE_SITES: {name} ({strand})',
                ),
                plot_components.Tracks(
                    tdata=output.splice_site_usage.filter_to_negative_strand(),
                    ylabel_template='SPLICE_SITE_USAGE: {biosample_name} ({strand})',
                ),
            ],
            interval=apol4_interval,
            title='Predicted SPLICE_SITE_USAGE for APOL4 gene in colon tissue',
        )

        # Save plot
        output_file = RESULTS_DIR / 'splice_site_usage.png'
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"✓ Plot saved: {output_file}")

        # Update results
        results['outputs']['splice_site_usage'] = {
            'file': str(output_file),
            'status': 'success'
        }
        results['status'] = 'success'

        # Save results metadata
        results_file = RESULTS_DIR / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)

        print(f"\n{'='*60}")
        print(f"✓ SPLICE_SITE_USAGE visualization completed successfully!")
        print(f"✓ Results saved to: {RESULTS_DIR}")
        print(f"✓ Metadata: {results_file}")
        print(f"{'='*60}")

        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        traceback.print_exc()

        # Save error results
        results['status'] = 'error'
        results['error'] = str(e)
        results_file = RESULTS_DIR / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)

        return 1


if __name__ == '__main__':
    sys.exit(main())
