#!/usr/bin/env python3
"""
AlphaGenome Analysis Workflow - TAL1 Locus Analysis
Analyzes non-coding T-ALL variants near the TAL1 locus
"""

import io
import itertools
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotnine as gg
from alphagenome.data import gene_annotation, genome, transcript as transcript_utils
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

# Configuration
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/analysis_workflow")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# API key from environment or .env file
API_KEY = os.environ.get("ALPHAGENOME_API_KEY")
if not API_KEY:
    # Try loading from .env file
    env_file = Path("/home/kyuwon/projects/alphagenome/.env")
    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                if line.startswith("ALPHAGENOME_API_KEY"):
                    API_KEY = line.split("=", 1)[1].strip().strip("'\"")
                    break
if not API_KEY:
    raise ValueError("ALPHAGENOME_API_KEY not found in environment or .env file")


def oncogenic_tal1_variants() -> pd.DataFrame:
    """Returns a dataframe of oncogenic T-ALL variants that affect TAL1."""
    variant_data = """
ID	CHROM	POS	REF	ALT	output	Study ID	Study Variant ID
Jurkat	chr1	47239296	C	CCGTTTCCTAACC	1	Mansour_2014
MOLT-3	chr1	47239296	C	ACC	1	Mansour_2014
Patient_1	chr1	47239296	C	AACG	1	Mansour_2014
Patient_2	chr1	47239291	CTAACC	TTTACCGTCTGTTAACGGC	1	Mansour_2014
Patient_3-5	chr1	47239296	C	ACG	1	Mansour_2014
Patient_6	chr1	47239296	C	ACC	1	Mansour_2014
Patient_7	chr1	47239295	AC	TCAAACTGGTAACC	1	Mansour_2014
Patient_8	chr1	47239296	C	AACC	1	Mansour_2014
new 3' enhancer 1	chr1	47212072	T	TGGGTAAACCGTCTGTTCAGCG	1	Smith_2023	UPNT802
new 3' enhancer 2	chr1	47212074	G	GAACGTT	1	Smith_2023	UPNT613
intergenic SNV 1	chr1	47230639	C	T	1	Liu_2020	SJALL043861_D1
intergenic SNV 2	chr1	47230639	C	T	1	Liu_2020	SJALL018373_D1
SJALL040467_D1	chr1	47239296	C	AACC	1	Liu_2020	SJALL040467_D1
PATBGC	chr1	47239296	C	AACC	1	Liu_2017	PATBGC
PATBTX	chr1	47239296	C	ACGGATATAACC	1	Liu_2017	PATBTX
PARJAY	chr1	47239296	C	ACGGAATTTCTAACC	1	Liu_2017	PARJAY
PARSJG	chr1	47239296	C	AACC	1	Liu_2017	PARSJG
PASYAJ	chr1	47239296	C	AACC	1	Liu_2017	PASYAJ
PATRAB	chr1	47239293	TTA	CTAACGG	1	Liu_2017	PATRAB
PAUBXP	chr1	47239296	C	ACC	1	Liu_2017	PAUBXP
PATENL	chr1	47239296	C	AACC	1	Liu_2017	PATENL
PARNXJ	chr1	47239296	C	ACG	1	Liu_2017	PARNXJ
PASXSI	chr1	47239296	C	AACC	1	Liu_2017	PASXSI
PASNEH	chr1	47239296	C	ACC	1	Liu_2017	PASNEH
PAUAFN	chr1	47239296	C	AACC	1	Liu_2017	PAUAFN
PARASZ	chr1	47239296	C	ACC	1	Liu_2017	PARASZ
PARWNW	chr1	47239296	C	ACC	1	Liu_2017	PARWNW
PASFKA	chr1	47239293	TTA	ACCGTTAATCAA	1	Liu_2017	PASFKA
PATEIT	chr1	47239296	C	AC	1	Liu_2017	PATEIT
PASMHF	chr1	47239296	C	AC	1	Liu_2017	PASMHF
PARJNX	chr1	47239296	C	AC	1	Liu_2017	PARJNX
PASYWF	chr1	47239296	C	AC	1	Liu_2017	PASYWF
"""
    return pd.read_table(io.StringIO(variant_data), sep='\t')


def generate_background_variants(variant: genome.Variant, max_number: int = 100) -> pd.DataFrame:
    """Generates background variants by random sequence generation."""
    nucleotides = np.array(list('ACGT'), dtype='<U1')

    def generate_unique_strings(n, max_number, random_seed=42):
        rng = np.random.default_rng(random_seed)
        if 4**n < max_number:
            raise ValueError('Cannot generate that many unique strings for the given length.')

        generated_strings = set()
        while len(generated_strings) < max_number:
            indices = rng.integers(0, 4, size=n)
            new_string = ''.join(nucleotides[indices])
            if new_string != variant.alternate_bases:
                generated_strings.add(new_string)
        return list(generated_strings)

    permutations = []
    if 4 ** len(variant.alternate_bases) < max_number:
        for p in itertools.product(nucleotides, repeat=len(variant.alternate_bases)):
            permutations.append(''.join(p))
    else:
        permutations = generate_unique_strings(len(variant.alternate_bases), max_number)

    ism_candidates = pd.DataFrame({
        'ID': ['mut_' + str(variant.position) + '_' + x for x in permutations],
        'CHROM': variant.chromosome,
        'POS': variant.position,
        'REF': variant.reference_bases,
        'ALT': permutations,
        'output': 0.0,
        'original_variant': variant.name,
    })
    return ism_candidates


def vcf_row_to_variant(vcf_row: pd.Series) -> genome.Variant:
    """Parse a row of a vcf df into a genome.Variant."""
    return genome.Variant(
        chromosome=str(vcf_row.CHROM),
        position=int(vcf_row.POS),
        reference_bases=vcf_row.REF,
        alternate_bases=vcf_row.ALT,
        name=vcf_row.ID,
    )


def inference_df(qtl_df: pd.DataFrame, input_sequence_length: int) -> pd.DataFrame:
    """Returns a pd.DataFrame with variants and intervals ready for inference."""
    df = []
    for _, row in qtl_df.iterrows():
        variant = vcf_row_to_variant(row)
        interval = genome.Interval(
            chromosome=row['CHROM'], start=row['POS'], end=row['POS']
        ).resize(input_sequence_length)

        df.append({
            'interval': interval,
            'variant': variant,
            'output': row['output'],
            'variant_id': row['ID'],
            'POS': row['POS'],
            'REF': row['REF'],
            'ALT': row['ALT'],
            'CHROM': row['CHROM'],
        })
    return pd.DataFrame(df)


def oncogenic_and_background_variants(input_sequence_length: int, number_of_background_variants: int = 20) -> pd.DataFrame:
    """Generates a dataframe of all variants for evaluation."""
    oncogenic_variants = oncogenic_tal1_variants()

    variants = []
    for vcf_row in oncogenic_variants.itertuples():
        variants.append(
            genome.Variant(
                chromosome=str(vcf_row.CHROM),
                position=int(vcf_row.POS),
                reference_bases=vcf_row.REF,
                alternate_bases=vcf_row.ALT,
                name=vcf_row.ID,
            )
        )

    background_variants = pd.concat([
        generate_background_variants(variant, number_of_background_variants)
        for variant in variants
    ])
    all_variants = pd.concat([oncogenic_variants, background_variants])
    return inference_df(all_variants, input_sequence_length=input_sequence_length)


def coarse_grained_mute_groups(eval_df):
    """Create variant groups for plotting."""
    grp = []
    for row in eval_df.itertuples():
        if row.POS >= 47239290:  # MUTE site
            if row.ALT_len > 4:
                grp.append('MUTE_other')
            else:
                grp.append('MUTE_' + str(row.ALT_len))
        else:
            grp.append(str(row.POS) + '_' + str(row.ALT_len))

    grp = pd.Series(grp)
    return pd.Categorical(grp, categories=sorted(grp.unique()), ordered=True)


def main():
    print("Starting AlphaGenome Analysis Workflow - TAL1 Locus")
    print(f"Results will be saved to: {RESULTS_DIR}")

    # Initialize model
    print("\n1. Initializing AlphaGenome DNA model...")
    dna_model = dna_client.create(API_KEY)

    # Load gene annotations
    print("\n2. Loading gene annotations from GENCODE...")
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    gtf_transcript = gene_annotation.filter_protein_coding(gtf)
    gtf_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcript)

    # Define TAL1 interval and ontology
    print("\n3. Setting up TAL1 analysis...")
    tal1_interval = genome.Interval(
        chromosome='chr1', start=47209255, end=47242023, strand='-'
    )
    ontology_terms = ['CL:0001059']  # CD34-positive common myeloid progenitor

    # Analyze single variant (Jurkat)
    print("\n4. Analyzing individual variant (Jurkat cell line)...")
    oncogenic_variants = oncogenic_tal1_variants()
    variant = vcf_row_to_variant(oncogenic_variants.iloc[0])
    print(f"   Variant: {variant}")

    print("   Making predictions (this may take a minute)...")
    output = dna_model.predict_variant(
        interval=tal1_interval.resize(2**20),
        variant=variant,
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.CHIP_HISTONE,
            dna_client.OutputType.DNASE,
        },
        ontology_terms=ontology_terms,
    )

    # Plot variant effect
    print("   Creating visualization...")
    transcripts = transcript_extractor.extract(tal1_interval)
    fig = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.Tracks(
                tdata=output.alternate.rna_seq.filter_to_nonpositive_strand()
                - output.reference.rna_seq.filter_to_nonpositive_strand(),
                ylabel_template='{biosample_name} ({strand})\n{name}',
                filled=True,
            ),
            plot_components.Tracks(
                tdata=output.alternate.dnase.filter_to_nonpositive_strand()
                - output.reference.dnase.filter_to_nonpositive_strand(),
                ylabel_template='{biosample_name} ({strand})\n{name}',
                filled=True,
            ),
            plot_components.Tracks(
                tdata=output.alternate.chip_histone.filter_to_nonpositive_strand()
                - output.reference.chip_histone.filter_to_nonpositive_strand(),
                ylabel_template='{biosample_name} ({strand})\n{name}',
                filled=True,
            ),
        ],
        annotations=[plot_components.VariantAnnotation([variant])],
        interval=tal1_interval,
        title=f'Effect of Jurkat variant on predicted outputs\n{variant}',
    )
    fig.set_size_inches(12, 8)
    fig.savefig(RESULTS_DIR / "jurkat_variant_effect.png", dpi=150, bbox_inches='tight')
    print(f"   Saved: {RESULTS_DIR / 'jurkat_variant_effect.png'}")

    # Compare oncogenic vs background variants
    print("\n5. Comparing oncogenic vs background variants...")
    print("   Generating variant sets (using 3 background variants per oncogenic variant)...")
    eval_df = oncogenic_and_background_variants(
        input_sequence_length=2**20,
        number_of_background_variants=3
    )

    eval_df['ALT_len'] = eval_df['ALT'].str.len()
    eval_df['variant_group'] = eval_df['POS'].astype(str) + '_' + eval_df['ALT_len'].astype(str)
    eval_df['output'] = eval_df['output'].fillna(0) != 0
    eval_df['coarse_grained_variant_group'] = coarse_grained_mute_groups(eval_df)

    print(f"   Total variants to score: {len(eval_df)}")
    print("   Scoring variants (this will take several minutes)...")
    scores = dna_model.score_variants(
        intervals=eval_df['interval'].to_list(),
        variants=eval_df['variant'].to_list(),
        variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
        max_workers=2,
    )

    # Extract TAL1 scores
    print("   Extracting TAL1 expression scores...")
    gene_index = scores[0][0].obs.query('gene_name == "TAL1"').index[0]
    cell_type_index = scores[0][0].var.query('ontology_curie == "CL:0001059"').index[0]

    def get_tal1_score_for_cd34_cells(score_data):
        return score_data[gene_index, cell_type_index].X[0, 0]

    eval_df['tal1_diff_in_cd34'] = [get_tal1_score_for_cd34_cells(x[0]) for x in scores]

    # Save results
    print("\n6. Saving results...")
    results = {
        'variant_scores': eval_df[[
            'variant_id', 'CHROM', 'POS', 'REF', 'ALT',
            'output', 'tal1_diff_in_cd34', 'variant_group'
        ]].to_dict('records'),
        'summary': {
            'total_variants': len(eval_df),
            'oncogenic_variants': int(eval_df['output'].sum()),
            'background_variants': int((~eval_df['output']).sum()),
            'mean_oncogenic_effect': float(eval_df[eval_df['output']]['tal1_diff_in_cd34'].mean()),
            'mean_background_effect': float(eval_df[~eval_df['output']]['tal1_diff_in_cd34'].mean()),
        }
    }

    with open(RESULTS_DIR / "variant_analysis_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    print(f"   Saved: {RESULTS_DIR / 'variant_analysis_results.json'}")

    # Create comparison plots
    print("\n7. Creating comparison plots...")
    plot_df = eval_df.loc[eval_df.REF != eval_df.ALT].copy()
    plot_df['variant'] = plot_df['variant'].astype(str)
    plot_df = plot_df[['variant', 'output', 'tal1_diff_in_cd34', 'coarse_grained_variant_group']].drop_duplicates()

    facet_title_by_group = {
        '47212072_22': 'chr1:47212072\n21 bp ins.',
        '47212074_7': 'chr1:47212074\n6 bp ins.',
        '47230639_1': 'chr1:47230639\nSNV',
        'MUTE_2': 'chr1:47239296\n1 bp ins.',
        'MUTE_3': 'chr1:47239296\n2 bp ins.',
        'MUTE_4': 'chr1:47239296\n3 bp ins.',
        'MUTE_other': 'chr1:47239296\n7-18 bp ins.',
    }

    for group in plot_df.coarse_grained_variant_group.unique():
        subplot_df = pd.concat([
            plot_df.assign(plot_group='density'),
            plot_df.assign(plot_group='rain')
        ])
        subplot_df = subplot_df[subplot_df.coarse_grained_variant_group == group]
        subplot_df = subplot_df[~((subplot_df.plot_group == 'density') & (subplot_df.output))]

        col_width = np.ptp(subplot_df.tal1_diff_in_cd34) / 200
        subplot_df['col_width'] = subplot_df['output'].map({True: 1.5 * col_width, False: 1.25 * col_width})

        plt_ = (
            gg.ggplot(subplot_df)
            + gg.aes(x='tal1_diff_in_cd34')
            + gg.geom_col(
                gg.aes(y=1, width='col_width', fill='output', x='tal1_diff_in_cd34', alpha='output'),
                data=subplot_df[subplot_df['plot_group'] == 'rain'],
            )
            + gg.geom_density(
                gg.aes(x='tal1_diff_in_cd34', fill='output'),
                data=subplot_df[subplot_df['plot_group'] == 'density'],
                color='white',
            )
            + gg.facet_wrap('~output + plot_group', nrow=1, scales='free_x')
            + gg.scale_alpha_manual({True: 1, False: 0.3})
            + gg.scale_fill_manual({True: '#FAA41A', False: 'gray'})
            + gg.labs(title=facet_title_by_group.get(group, group))
            + gg.theme_minimal()
            + gg.geom_vline(xintercept=0, linetype='dotted')
            + gg.theme(
                figure_size=(4, 3),
                legend_position='none',
                axis_text_x=gg.element_blank(),
                panel_grid_major_x=gg.element_blank(),
                panel_grid_minor_x=gg.element_blank(),
                strip_text=gg.element_blank(),
                axis_title_y=gg.element_blank(),
                axis_title_x=gg.element_blank(),
                plot_title=gg.element_text(size=10),
            )
            + gg.scale_y_reverse()
            + gg.coord_flip()
        )

        filename = f"comparison_{group}.png"
        plt_.save(RESULTS_DIR / filename, dpi=150)
        print(f"   Saved: {RESULTS_DIR / filename}")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nResults directory: {RESULTS_DIR}")
    print("\nGenerated files:")
    print("  - jurkat_variant_effect.png: Individual variant effect visualization")
    print("  - variant_analysis_results.json: Quantitative results")
    print("  - comparison_*.png: Oncogenic vs background variant comparisons")
    print("\nSummary:")
    print(f"  Total variants analyzed: {results['summary']['total_variants']}")
    print(f"  Oncogenic variants: {results['summary']['oncogenic_variants']}")
    print(f"  Background variants: {results['summary']['background_variants']}")
    print(f"  Mean oncogenic TAL1 effect: {results['summary']['mean_oncogenic_effect']:.4f}")
    print(f"  Mean background TAL1 effect: {results['summary']['mean_background_effect']:.4f}")


if __name__ == "__main__":
    main()
