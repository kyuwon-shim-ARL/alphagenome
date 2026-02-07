#!/usr/bin/env python3
"""
Variant Scoring CLI Script
CLI-based variant scoring and visualization, replacing Colab-dependent variant_scoring_ui.ipynb

Usage:
    python scripts/run_variant_scoring_cli.py \
        --chr chr22 --pos 36201698 --ref A --alt C \
        --outputs rna_seq,cage,dnase \
        --ontology UBERON:0001157

For full options:
    python scripts/run_variant_scoring_cli.py --help
"""

import argparse
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
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

load_dotenv()

# Constants
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/variant_scoring_cli")
API_DELAY = 1.0

# Output type mapping
OUTPUT_TYPE_MAP = {
    'rna_seq': dna_client.OutputType.RNA_SEQ,
    'cage': dna_client.OutputType.CAGE,
    'dnase': dna_client.OutputType.DNASE,
    'atac': dna_client.OutputType.ATAC,
    'chip_histone': dna_client.OutputType.CHIP_HISTONE,
    'chip_tf': dna_client.OutputType.CHIP_TF,
    'splice_sites': dna_client.OutputType.SPLICE_SITES,
    'splice_site_usage': dna_client.OutputType.SPLICE_SITE_USAGE,
    'splice_junctions': dna_client.OutputType.SPLICE_JUNCTIONS,
    'contact_maps': dna_client.OutputType.CONTACT_MAPS,
    'procap': dna_client.OutputType.PROCAP,
}

# Gene annotation URLs
HG38_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='CLI-based variant scoring and visualization using AlphaGenome',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Score a variant with default outputs
    python scripts/run_variant_scoring_cli.py --chr chr22 --pos 36201698 --ref A --alt C

    # Score with specific outputs and ontology terms
    python scripts/run_variant_scoring_cli.py \\
        --chr chr22 --pos 36201698 --ref A --alt C \\
        --outputs rna_seq,cage,dnase,atac,chip_histone \\
        --ontology UBERON:0001157,EFO:0001187

    # Use 500KB sequence length and custom plot interval
    python scripts/run_variant_scoring_cli.py \\
        --chr chr22 --pos 36201698 --ref A --alt C \\
        --sequence-length 500KB --interval-width 32768
        """
    )

    # Required arguments
    parser.add_argument('--chr', required=True, help='Chromosome (e.g., chr22)')
    parser.add_argument('--pos', type=int, required=True, help='Position (e.g., 36201698)')
    parser.add_argument('--ref', required=True, help='Reference base(s) (e.g., A)')
    parser.add_argument('--alt', required=True, help='Alternate base(s) (e.g., C)')

    # Optional arguments
    parser.add_argument('--outputs', default='rna_seq,dnase',
                        help='Comma-separated output types (default: rna_seq,dnase). '
                             f'Available: {",".join(OUTPUT_TYPE_MAP.keys())}')
    parser.add_argument('--ontology', default=None,
                        help='Comma-separated ontology terms (e.g., UBERON:0001157,EFO:0001187). '
                             'If not specified, uses all available.')
    parser.add_argument('--sequence-length', default='1MB',
                        choices=['2KB', '16KB', '100KB', '500KB', '1MB'],
                        help='Sequence length around variant (default: 1MB)')
    parser.add_argument('--interval-width', type=int, default=32768,
                        help='Visualization interval width in bp (default: 32768)')
    parser.add_argument('--no-visualize', action='store_true',
                        help='Skip visualization, only generate scores')
    parser.add_argument('--filter-strand', choices=['+', '-'],
                        help='Filter tracks to specific DNA strand')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Custom output directory (default: results/variant_scoring_cli)')

    return parser.parse_args()


def get_api_key():
    """Retrieve API key from environment"""
    api_key = os.getenv('ALPHAGENOME_API_KEY')
    if not api_key:
        raise ValueError(
            "ALPHAGENOME_API_KEY environment variable not set. "
            "Please set it before running this script."
        )
    return api_key


def score_variant(dna_model, variant, interval, organism):
    """Score variant using all recommended scorers"""
    print("\n[Step 1] Scoring variant...")

    variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
        organism=organism,
    )

    df_scores = variant_scorers.tidy_scores(variant_scores)
    print(f"  Generated {len(df_scores)} score rows")

    return df_scores


def visualize_variant(dna_model, variant, interval, output_types, ontology_terms,
                      filter_strand, interval_width, results_dir, transcript_extractor):
    """Generate visualization plots for variant effects"""
    print("\n[Step 2] Generating visualizations...")

    # Predict variant effects
    print(f"  Requesting outputs: {[ot.name for ot in output_types]}")
    if ontology_terms:
        print(f"  Ontology terms: {ontology_terms}")

    output = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=output_types,
        ontology_terms=ontology_terms if ontology_terms else None,
    )
    time.sleep(API_DELAY)

    ref, alt = output.reference, output.alternate

    # Filter by strand if requested
    if filter_strand == '+':
        ref = ref.filter_to_strand(strand='+')
        alt = alt.filter_to_strand(strand='+')
    elif filter_strand == '-':
        ref = ref.filter_to_strand(strand='-')
        alt = alt.filter_to_strand(strand='-')

    # Get transcripts for annotation
    transcripts = transcript_extractor.extract(interval)

    # Plot interval centered on variant
    plot_interval = variant.reference_interval.resize(interval_width)

    # Color scheme
    ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}

    # Generate plots for each output type
    plots_generated = []

    output_map = {
        dna_client.OutputType.RNA_SEQ: ('rna_seq', ref.rna_seq, alt.rna_seq, 'RNA_SEQ: {biosample_name} ({strand})'),
        dna_client.OutputType.CAGE: ('cage', ref.cage, alt.cage, 'CAGE: {biosample_name} ({strand})'),
        dna_client.OutputType.DNASE: ('dnase', ref.dnase, alt.dnase, 'DNASE: {biosample_name} ({strand})'),
        dna_client.OutputType.ATAC: ('atac', ref.atac, alt.atac, 'ATAC: {biosample_name} ({strand})'),
        dna_client.OutputType.CHIP_HISTONE: ('chip_histone', ref.chip_histone, alt.chip_histone,
                                              'CHIP_HISTONE: {biosample_name} ({strand})\n{histone_mark}'),
        dna_client.OutputType.CHIP_TF: ('chip_tf', ref.chip_tf, alt.chip_tf,
                                         'CHIP_TF: {biosample_name} ({strand})\n{transcription_factor}'),
        dna_client.OutputType.SPLICE_SITES: ('splice_sites', ref.splice_sites, alt.splice_sites,
                                              'SPLICE_SITES: {name} ({strand})'),
    }

    for output_type in output_types:
        if output_type not in output_map:
            continue

        name, ref_data, alt_data, ylabel_template = output_map[output_type]

        if ref_data is None or alt_data is None:
            print(f"  Skipping {name}: no data available")
            continue

        if ref_data.values.shape[-1] == 0:
            print(f"  Skipping {name}: no tracks for specified ontology/strand")
            continue

        print(f"  Generating plot for {name}...")

        components = [
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.OverlaidTracks(
                tdata={'REF': ref_data, 'ALT': alt_data},
                colors=ref_alt_colors,
                ylabel_template=ylabel_template,
            ),
        ]

        fig = plot_components.plot(
            components=components,
            interval=plot_interval,
            annotations=[
                plot_components.VariantAnnotation([variant]),
            ],
            title=f'Variant Effect: {variant} - {name.upper()}',
        )

        filepath = results_dir / f'plot_{name}.png'
        fig.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close(fig)
        plots_generated.append(str(filepath))
        print(f"    Saved: {filepath}")

    return plots_generated


def main():
    args = parse_args()

    # Setup results directory
    results_dir = Path(args.output_dir) if args.output_dir else RESULTS_DIR
    results_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Initialize
        print("=" * 70)
        print("AlphaGenome Variant Scoring CLI")
        print("=" * 70)

        api_key = get_api_key()
        dna_model = dna_client.create(api_key)

        # Create variant
        variant = genome.Variant(
            chromosome=args.chr,
            position=args.pos,
            reference_bases=args.ref,
            alternate_bases=args.alt,
        )
        print(f"\nVariant: {variant}")

        # Set sequence length
        seq_length_key = f'SEQUENCE_LENGTH_{args.sequence_length}'
        sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[seq_length_key]
        interval = variant.reference_interval.resize(sequence_length)
        print(f"Interval: {interval}")

        # Parse output types
        output_type_names = [o.strip().lower() for o in args.outputs.split(',')]
        output_types = []
        for name in output_type_names:
            if name in OUTPUT_TYPE_MAP:
                output_types.append(OUTPUT_TYPE_MAP[name])
            else:
                print(f"Warning: Unknown output type '{name}', skipping")
        print(f"Output types: {[ot.name for ot in output_types]}")

        # Parse ontology terms
        ontology_terms = None
        if args.ontology:
            ontology_terms = [o.strip() for o in args.ontology.split(',')]
            print(f"Ontology terms: {ontology_terms}")

        # Score variant
        df_scores = score_variant(dna_model, variant, interval, dna_client.Organism.HOMO_SAPIENS)

        # Save scores
        scores_path = results_dir / 'variant_scores.csv'
        df_scores.to_csv(scores_path, index=False)
        print(f"\n  Saved scores to: {scores_path}")

        # Create summary
        summary = {
            'variant': str(variant),
            'chromosome': args.chr,
            'position': args.pos,
            'reference': args.ref,
            'alternate': args.alt,
            'interval': str(interval),
            'sequence_length': args.sequence_length,
            'total_scores': len(df_scores),
            'score_statistics': {
                'mean_raw_score': float(df_scores['raw_score'].mean()),
                'median_raw_score': float(df_scores['raw_score'].median()),
                'std_raw_score': float(df_scores['raw_score'].std()),
                'min_raw_score': float(df_scores['raw_score'].min()),
                'max_raw_score': float(df_scores['raw_score'].max()),
            },
            'unique_scorers': df_scores['variant_scorer'].nunique() if 'variant_scorer' in df_scores.columns else 0,
        }

        # Optional visualization
        plots_generated = []
        if not args.no_visualize and output_types:
            # Load gene annotations
            print("\n  Loading gene annotations...")
            gtf = pd.read_feather(HG38_GTF_FEATHER)
            gtf_transcript = gene_annotation.filter_transcript_support_level(
                gene_annotation.filter_protein_coding(gtf), ['1']
            )
            transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)

            plots_generated = visualize_variant(
                dna_model, variant, interval, output_types, ontology_terms,
                args.filter_strand, args.interval_width, results_dir, transcript_extractor
            )
            summary['plots'] = plots_generated

        # Save summary
        summary_path = results_dir / 'variant_summary.json'
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"  Saved summary to: {summary_path}")

        # Print results summary
        print("\n" + "=" * 70)
        print("VARIANT SCORING COMPLETE")
        print("=" * 70)
        print(f"\nVariant: {variant}")
        print(f"Total scores: {len(df_scores)}")
        print(f"Plots generated: {len(plots_generated)}")
        print(f"\nResults directory: {results_dir}")
        print(f"  - variant_scores.csv ({len(df_scores)} rows)")
        print(f"  - variant_summary.json")
        for plot in plots_generated:
            print(f"  - {Path(plot).name}")

        # Show top scores
        print("\nTop 5 positive effect scores:")
        top_cols = ['variant_scorer', 'ontology_curie', 'raw_score', 'quantile_score']
        top_cols = [c for c in top_cols if c in df_scores.columns]
        print(df_scores.nlargest(5, 'raw_score')[top_cols].to_string(index=False))

        print("\nTop 5 negative effect scores:")
        print(df_scores.nsmallest(5, 'raw_score')[top_cols].to_string(index=False))

        return 0

    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
