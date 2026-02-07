#!/usr/bin/env python3
"""
CHIP_TF Analysis Script
Generates transcription factor binding profile visualizations for K562 and HepG2 cell lines.

This script complements the visualization_modality_tour by adding CHIP_TF analysis
which was skipped in the original tutorial execution.

Output:
    results/chip_tf_analysis/
    ├── chip_tf_k562_ctcf.png         # K562 CTCF binding profile
    ├── chip_tf_hepg2_ctcf.png        # HepG2 CTCF binding profile
    ├── chip_tf_ctcf_rad21_coloc.png  # CTCF-RAD21 co-localization comparison
    ├── chip_tf_multi_tf.png          # Multiple TF comparison
    └── results.json                   # Metadata and results summary
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
import numpy as np
from dotenv import load_dotenv

from alphagenome.data import gene_annotation, genome, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

load_dotenv()

# Configuration
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/chip_tf_analysis")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

API_DELAY = 2.0

# Cell line ontology codes
K562_ONTOLOGY = "EFO:0002067"    # K562 chronic myelogenous leukemia cell line
HEPG2_ONTOLOGY = "EFO:0001187"  # HepG2 hepatocellular carcinoma cell line

# Transcription factors to analyze
TARGET_TFS = ["CTCF", "RAD21", "POLR2A", "EP300"]

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


def filter_to_tfs(chip_tf_data, tf_names):
    """Filter CHIP_TF TrackData to specific transcription factors"""
    if chip_tf_data is None:
        return None

    metadata = chip_tf_data.metadata
    if 'transcription_factor' not in metadata.columns:
        print(f"  Warning: 'transcription_factor' column not found in metadata")
        return chip_tf_data

    # Case-insensitive matching
    tf_mask = metadata['transcription_factor'].str.upper().isin([tf.upper() for tf in tf_names])

    if not tf_mask.any():
        print(f"  Warning: No tracks found for TFs: {tf_names}")
        return None

    filtered_indices = metadata.index[tf_mask]
    return chip_tf_data.select_tracks_by_index(filtered_indices)


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
        print("AlphaGenome CHIP_TF Analysis")
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
        results['metadata']['cell_lines'] = {
            'K562': K562_ONTOLOGY,
            'HepG2': HEPG2_ONTOLOGY,
        }
        results['metadata']['target_tfs'] = TARGET_TFS

        # ==================== 1. K562 CTCF Binding Profile ====================
        print("\n[1/4] Generating K562 CTCF binding profile...")

        output_k562 = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.CHIP_TF},
            ontology_terms=[K562_ONTOLOGY],
        )
        time.sleep(API_DELAY)

        k562_ctcf = filter_to_tfs(output_k562.chip_tf, ['CTCF'])

        if k562_ctcf is not None and k562_ctcf.values.shape[-1] > 0:
            fig = plot_components.plot(
                [
                    plot_components.TranscriptAnnotation(apol_transcripts),
                    plot_components.Tracks(
                        tdata=k562_ctcf,
                        ylabel_template='CHIP_TF: K562\n{transcription_factor} ({strand})',
                        filled=True,
                    ),
                ],
                interval=interval.resize(200_000),  # Focus on 200kb
                title='CTCF Binding Profile in K562 (Leukemia Cell Line)',
            )
            save_plot(fig, 'chip_tf_k562_ctcf', results['plots'])
        else:
            print("  Warning: No CTCF data available for K562")
            results['plots']['chip_tf_k562_ctcf'] = {'status': 'no_data'}

        # ==================== 2. HepG2 CTCF Binding Profile ====================
        print("\n[2/4] Generating HepG2 CTCF binding profile...")

        output_hepg2 = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.CHIP_TF},
            ontology_terms=[HEPG2_ONTOLOGY],
        )
        time.sleep(API_DELAY)

        hepg2_ctcf = filter_to_tfs(output_hepg2.chip_tf, ['CTCF'])

        if hepg2_ctcf is not None and hepg2_ctcf.values.shape[-1] > 0:
            fig = plot_components.plot(
                [
                    plot_components.TranscriptAnnotation(apol_transcripts),
                    plot_components.Tracks(
                        tdata=hepg2_ctcf,
                        ylabel_template='CHIP_TF: HepG2\n{transcription_factor} ({strand})',
                        filled=True,
                    ),
                ],
                interval=interval.resize(200_000),
                title='CTCF Binding Profile in HepG2 (Liver Cancer Cell Line)',
            )
            save_plot(fig, 'chip_tf_hepg2_ctcf', results['plots'])
        else:
            print("  Warning: No CTCF data available for HepG2")
            results['plots']['chip_tf_hepg2_ctcf'] = {'status': 'no_data'}

        # ==================== 3. CTCF-RAD21 Co-localization ====================
        print("\n[3/4] Generating CTCF-RAD21 co-localization comparison...")

        # Use K562 for co-localization analysis
        k562_ctcf_rad21 = filter_to_tfs(output_k562.chip_tf, ['CTCF', 'RAD21'])

        if k562_ctcf_rad21 is not None and k562_ctcf_rad21.values.shape[-1] > 0:
            # Define TF colors
            tf_colors = []
            for tf in k562_ctcf_rad21.metadata['transcription_factor']:
                if 'CTCF' in tf.upper():
                    tf_colors.append('#e41a1c')  # Red for CTCF
                elif 'RAD21' in tf.upper():
                    tf_colors.append('#377eb8')  # Blue for RAD21
                else:
                    tf_colors.append('#999999')

            fig = plot_components.plot(
                [
                    plot_components.TranscriptAnnotation(apol_transcripts),
                    plot_components.Tracks(
                        tdata=k562_ctcf_rad21,
                        ylabel_template='CHIP_TF: K562\n{transcription_factor} ({strand})',
                        filled=True,
                        track_colors=tf_colors,
                    ),
                ],
                interval=interval.resize(200_000),
                title='CTCF-RAD21 Co-localization in K562\n(Red: CTCF, Blue: RAD21)',
            )
            save_plot(fig, 'chip_tf_ctcf_rad21_coloc', results['plots'])
        else:
            print("  Warning: No CTCF/RAD21 data available for co-localization")
            results['plots']['chip_tf_ctcf_rad21_coloc'] = {'status': 'no_data'}

        # ==================== 4. Multi-TF Comparison ====================
        print("\n[4/4] Generating multi-TF comparison...")

        k562_multi_tf = filter_to_tfs(output_k562.chip_tf, TARGET_TFS)

        if k562_multi_tf is not None and k562_multi_tf.values.shape[-1] > 0:
            # Color by TF
            tf_color_map = {
                'CTCF': '#e41a1c',
                'RAD21': '#377eb8',
                'POLR2A': '#4daf4a',
                'EP300': '#984ea3',
            }
            multi_tf_colors = []
            for tf in k562_multi_tf.metadata['transcription_factor']:
                tf_upper = tf.upper()
                color = '#999999'
                for key, val in tf_color_map.items():
                    if key in tf_upper:
                        color = val
                        break
                multi_tf_colors.append(color)

            fig = plot_components.plot(
                [
                    plot_components.TranscriptAnnotation(apol_transcripts),
                    plot_components.Tracks(
                        tdata=k562_multi_tf,
                        ylabel_template='CHIP_TF: K562\n{transcription_factor} ({strand})',
                        filled=True,
                        track_colors=multi_tf_colors,
                    ),
                ],
                interval=interval.resize(200_000),
                title='Multi-TF Binding Profiles in K562\n(CTCF, RAD21, POLR2A, EP300)',
            )
            save_plot(fig, 'chip_tf_multi_tf', results['plots'])

            # Record available TFs
            available_tfs = k562_multi_tf.metadata['transcription_factor'].unique().tolist()
            results['metadata']['available_tfs_k562'] = available_tfs
            print(f"  Available TFs in K562: {available_tfs}")
        else:
            print("  Warning: No multi-TF data available")
            results['plots']['chip_tf_multi_tf'] = {'status': 'no_data'}

        # Save results summary
        results['summary'] = {
            'total_plots': len([p for p in results['plots'].values() if p.get('status') == 'success']),
            'cell_lines_analyzed': ['K562', 'HepG2'],
            'target_tfs': TARGET_TFS,
        }

        results_file = RESULTS_DIR / 'results.json'
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        # Print summary
        print("\n" + "=" * 70)
        print("CHIP_TF ANALYSIS COMPLETE")
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
