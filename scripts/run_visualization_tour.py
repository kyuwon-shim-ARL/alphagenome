#!/usr/bin/env python3
"""
AlphaGenome Visualization Modality Tour
Executes key visualizations from the tutorial notebook
"""

import json
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import gene_annotation, genome, track_data
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

# Load environment variables from .env file
load_dotenv()

# Configuration
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/visualization_tour")
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

def save_plot(fig, name, results):
    """Save plot and record in results"""
    filepath = RESULTS_DIR / f"{name}.png"
    fig.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    results[name] = {
        'file': str(filepath),
        'status': 'success'
    }
    print(f"✓ Saved: {filepath}")

def main():
    results = {}

    try:
        print("Initializing AlphaGenome client...")
        api_key = get_api_key()
        dna_model = dna_client.create(api_key)

        # Load output metadata
        print("Loading output metadata...")
        output_metadata = dna_model.output_metadata(
            organism=dna_client.Organism.HOMO_SAPIENS
        )

        # Load gene annotations
        print("Loading gene annotations...")
        gtf = pd.read_feather(
            'https://storage.googleapis.com/alphagenome/reference/gencode/'
            'hg38/gencode.v46.annotation.gtf.gz.feather'
        )

        gtf_transcript = gene_annotation.filter_transcript_support_level(
            gene_annotation.filter_protein_coding(gtf), ['1']
        )

        gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(
            gtf_transcript
        )

        # Define genomic interval (chr22 around APOL4 gene)
        interval = genome.Interval('chr22', 36_150_498, 36_252_898).resize(
            dna_client.SEQUENCE_LENGTH_1MB
        )

        print(f"Working with interval: {interval}")

        # ==================== 1. RNA_SEQ Visualization ====================
        print("\n[1/7] Generating RNA_SEQ visualization...")
        ontology_terms = [
            'UBERON:0001159',  # Colon - Sigmoid
            'UBERON:0001155',  # Colon - Transverse
        ]

        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.RNA_SEQ},
            ontology_terms=ontology_terms,
        )
        time.sleep(API_DELAY)

        longest_transcripts = gene_annotation.filter_to_longest_transcript(
            gtf_transcript[
                (gtf_transcript['Chromosome'] == 'chr22') &
                (gtf_transcript['Start'] >= interval.start) &
                (gtf_transcript['End'] <= interval.end)
            ]
        )

        from alphagenome.data import transcript
        transcript_extractor = transcript.TranscriptExtractor(
            gtf_transcript
        )
        longest_transcripts_obj = transcript_extractor.extract(interval)
        longest_transcripts_obj = [
            t for t in longest_transcripts_obj
            if t.info.get('gene_name') in ['APOL1', 'APOL2', 'APOL3', 'APOL4', 'APOL5', 'APOL6']
        ]

        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(longest_transcripts_obj),
                plot_components.Tracks(
                    tdata=output.rna_seq,
                    ylabel_template='RNA_SEQ: {biosample_name} ({strand})',
                ),
            ],
            interval=interval,
            title='Predicted RNA Expression (RNA_SEQ) for colon tissue',
        )
        save_plot(fig, '01_rna_seq', results)

        # ==================== 2. CAGE Visualization ====================
        print("\n[2/7] Generating CAGE visualization...")
        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.CAGE},
            ontology_terms=ontology_terms,
        )
        time.sleep(API_DELAY)

        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(longest_transcripts_obj),
                plot_components.Tracks(
                    tdata=output.cage,
                    ylabel_template='CAGE: {biosample_name} ({strand})',
                ),
            ],
            interval=interval,
            title='Predicted CAGE (transcription start sites) for colon tissue',
        )
        save_plot(fig, '02_cage', results)

        # ==================== 3. DNASE Visualization ====================
        print("\n[3/7] Generating DNASE visualization...")
        ontology_terms_intestine = [
            'UBERON:0001155',  # Colon - Transverse
            'UBERON:0001159',  # Colon - Sigmoid
        ]

        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.DNASE},
            ontology_terms=ontology_terms_intestine,
        )
        time.sleep(API_DELAY)

        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(longest_transcripts_obj),
                plot_components.Tracks(
                    tdata=output.dnase,
                    ylabel_template='DNASE: {biosample_name} ({strand})',
                ),
            ],
            interval=interval,
            title='Predicted chromatin accessibility (DNASE) for colon tissue',
        )
        save_plot(fig, '03_dnase', results)

        # ==================== 4. ATAC Visualization ====================
        print("\n[4/7] Generating ATAC visualization...")
        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.ATAC},
            ontology_terms=ontology_terms_intestine,
        )
        time.sleep(API_DELAY)

        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(longest_transcripts_obj),
                plot_components.Tracks(
                    tdata=output.atac,
                    ylabel_template='ATAC: {biosample_name} ({strand})',
                ),
            ],
            interval=interval,
            title='Predicted chromatin accessibility (ATAC) for colon tissue',
        )
        save_plot(fig, '04_atac', results)

        # ==================== 5. CHIP_HISTONE Visualization ====================
        print("\n[5/7] Generating CHIP_HISTONE visualization...")
        ontology_terms_colon = [
            'UBERON:0000317',
            'UBERON:0001157',
            'UBERON:0001159',
        ]

        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.CHIP_HISTONE},
            ontology_terms=ontology_terms_colon,
        )
        time.sleep(API_DELAY)

        # Reorder by histone mark
        reordered_chip_histone = output.chip_histone.select_tracks_by_index(
            output.chip_histone.metadata.sort_values('histone_mark').index
        )

        histone_to_color = {
            'H3K27AC': '#e41a1c',
            'H3K36ME3': '#ff7f00',
            'H3K4ME1': '#377eb8',
            'H3K4ME3': '#984ea3',
            'H3K9AC': '#4daf4a',
            'H3K27ME3': '#ffc0cb',
        }

        track_colors = (
            reordered_chip_histone.metadata['histone_mark']
            .map(lambda x: histone_to_color.get(x.upper(), '#000000'))
            .values
        )

        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(longest_transcripts_obj),
                plot_components.Tracks(
                    tdata=reordered_chip_histone,
                    ylabel_template='CHIP_HISTONE: {biosample_name}\n{histone_mark}',
                    filled=True,
                    track_colors=track_colors,
                ),
            ],
            interval=interval,
            title='Predicted histone modification markers in colon tissue',
            despine_keep_bottom=True,
        )
        save_plot(fig, '05_chip_histone', results)

        # ==================== 6. SPLICE_SITES and SPLICE_JUNCTIONS ====================
        print("\n[6/7] Generating SPLICE_SITES and SPLICE_JUNCTIONS visualization...")
        ontology_terms_splice = [
            'UBERON:0001157',
            'UBERON:0001159',
        ]

        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={
                dna_client.OutputType.SPLICE_SITES,
                dna_client.OutputType.SPLICE_JUNCTIONS,
            },
            ontology_terms=ontology_terms_splice,
        )
        time.sleep(API_DELAY)

        # Focus on APOL4 interval
        apol4_interval = gene_annotation.get_gene_interval(gtf, gene_symbol='APOL4')
        apol4_interval.resize_inplace(apol4_interval.width + 1000)

        all_transcripts = transcript_extractor.extract(interval)

        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(all_transcripts),
                plot_components.Tracks(
                    tdata=output.splice_sites.filter_to_negative_strand(),
                    ylabel_template='SPLICE_SITES: {name} ({strand})',
                ),
                plot_components.Sashimi(
                    output.splice_junctions
                    .filter_to_strand('-')
                    .filter_by_tissue('Colon_Transverse'),
                    ylabel_template='SPLICE_JUNCTIONS: {biosample_name} ({strand})',
                ),
            ],
            interval=apol4_interval,
            title='Predicted splicing effects for colon tissue',
        )
        save_plot(fig, '06_splice', results)

        # ==================== 7. CONTACT_MAPS Visualization ====================
        print("\n[7/7] Generating CONTACT_MAPS visualization...")
        ontology_terms_contact = [
            'EFO:0002824',  # HCT116 colon carcinoma cell line
        ]

        output = dna_model.predict_interval(
            interval=interval,
            requested_outputs={dna_client.OutputType.CONTACT_MAPS},
            ontology_terms=ontology_terms_contact,
        )
        time.sleep(API_DELAY)

        fig = plot_components.plot(
            [
                plot_components.TranscriptAnnotation(longest_transcripts_obj),
                plot_components.ContactMaps(
                    tdata=output.contact_maps,
                    ylabel_template='{biosample_name}\n{name}',
                    cmap='autumn_r',
                    vmax=1.0,
                ),
            ],
            interval=interval,
            title='Predicted contact maps (3D chromatin interaction)',
        )
        save_plot(fig, '07_contact_maps', results)

        # Save results summary
        results_file = RESULTS_DIR / 'results.json'
        with open(results_file, 'w') as f:
            json.dump({
                'interval': str(interval),
                'visualizations': results,
                'summary': {
                    'total': len(results),
                    'successful': sum(1 for r in results.values() if r['status'] == 'success'),
                }
            }, f, indent=2)

        print(f"\n{'='*60}")
        print(f"✓ All visualizations completed successfully!")
        print(f"✓ Results saved to: {RESULTS_DIR}")
        print(f"✓ Summary: {results_file}")
        print(f"{'='*60}")

        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

        # Save partial results
        if results:
            results_file = RESULTS_DIR / 'results.json'
            with open(results_file, 'w') as f:
                json.dump({
                    'error': str(e),
                    'partial_results': results
                }, f, indent=2)

        return 1

if __name__ == '__main__':
    sys.exit(main())
