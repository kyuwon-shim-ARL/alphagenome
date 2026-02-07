#!/usr/bin/env python3
"""
AlphaGenome Tissue Ontology Mapping Tutorial

This script explores tissue and cell type ontologies used in AlphaGenome predictions.
It demonstrates:
- Navigating UBERON (anatomy) ontology
- Exploring CL (cell type) ontology
- Understanding EFO (experimental conditions) ontology
- Searching ontology terms for tissues/cells
- Listing available ontologies for predictions
"""

import os
import json
import pandas as pd
from pathlib import Path
from alphagenome.models import dna_client

# Configuration
RESULTS_DIR = Path("/home/kyuwon/projects/alphagenome/results/tissue_ontology")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# API Key (from environment)
API_KEY = os.getenv("ALPHAGENOME_API_KEY")
if not API_KEY:
    raise ValueError("ALPHAGENOME_API_KEY environment variable not set")


def explore_output_metadata(dna_model):
    """
    Explore available output metadata for human and mouse.

    Returns:
        dict: Contains metadata DataFrames and track counts
    """
    print("=" * 80)
    print("EXPLORING OUTPUT METADATA")
    print("=" * 80)

    # Get human metadata
    print("\nFetching human output metadata...")
    human_metadata = dna_model.output_metadata(
        dna_client.Organism.HOMO_SAPIENS
    ).concatenate()

    print(f"Total human tracks: {len(human_metadata)}")
    print(f"Columns: {list(human_metadata.columns)}")
    print(f"\nFirst few entries:")
    print(human_metadata.head(10))

    # Get mouse metadata
    print("\n" + "-" * 80)
    print("Fetching mouse output metadata...")
    mouse_metadata = dna_model.output_metadata(
        dna_client.Organism.MUS_MUSCULUS
    ).concatenate()

    print(f"Total mouse tracks: {len(mouse_metadata)}")

    # Count tracks per output type
    print("\n" + "-" * 80)
    print("TRACKS PER OUTPUT TYPE")
    print("-" * 80)

    human_tracks = (
        human_metadata
        .groupby('output_type')
        .size()
        .rename('# Human tracks')
    )

    mouse_tracks = (
        mouse_metadata
        .groupby('output_type')
        .size()
        .rename('# Mouse tracks')
    )

    track_counts = pd.concat([human_tracks, mouse_tracks], axis=1).astype(pd.Int64Dtype())
    print(track_counts)

    return {
        'human_metadata': human_metadata,
        'mouse_metadata': mouse_metadata,
        'track_counts': track_counts
    }


def extract_ontology_terms(metadata_df):
    """
    Extract unique ontology terms from metadata.

    Args:
        metadata_df: DataFrame with ontology information

    Returns:
        dict: Ontology types and their terms
    """
    print("\n" + "=" * 80)
    print("EXTRACTING ONTOLOGY TERMS")
    print("=" * 80)

    ontologies = {}

    # Extract from ontology_curie column if exists
    if 'ontology_curie' in metadata_df.columns:
        terms = metadata_df['ontology_curie'].dropna().unique()

        # Group by ontology prefix
        for term in terms:
            if ':' in term:
                prefix = term.split(':')[0]
                if prefix not in ontologies:
                    ontologies[prefix] = []
                ontologies[prefix].append(term)

    # Print summary
    for ont_type, terms in sorted(ontologies.items()):
        print(f"\n{ont_type}: {len(terms)} unique terms")
        print(f"  Examples: {terms[:5]}")

    return ontologies


def search_tissue_terms(metadata_df, search_terms):
    """
    Search for specific tissue/cell types in metadata.

    Args:
        metadata_df: Metadata DataFrame
        search_terms: List of terms to search for

    Returns:
        dict: Search results
    """
    print("\n" + "=" * 80)
    print("SEARCHING TISSUE/CELL TERMS")
    print("=" * 80)

    results = {}

    for term in search_terms:
        print(f"\nSearching for: '{term}'")

        # Search in all string columns
        mask = metadata_df.apply(
            lambda col: col.astype(str).str.contains(term, case=False, na=False)
        ).any(axis=1)

        matches = metadata_df[mask]

        print(f"  Found {len(matches)} matches")

        if len(matches) > 0:
            # Show unique ontology terms
            if 'ontology_curie' in matches.columns:
                unique_terms = matches['ontology_curie'].dropna().unique()
                print(f"  Unique ontology terms: {list(unique_terms[:5])}")

            # Show output types
            if 'output_type' in matches.columns:
                output_types = matches['output_type'].value_counts()
                print(f"  Output types:")
                for otype, count in output_types.items():
                    print(f"    {otype}: {count}")

        results[term] = {
            'count': len(matches),
            'matches': matches
        }

    return results


def analyze_ontology_coverage(metadata_df):
    """
    Analyze which ontologies are used for different output types.

    Args:
        metadata_df: Metadata DataFrame

    Returns:
        DataFrame: Coverage matrix
    """
    print("\n" + "=" * 80)
    print("ONTOLOGY COVERAGE BY OUTPUT TYPE")
    print("=" * 80)

    if 'ontology_curie' not in metadata_df.columns or 'output_type' not in metadata_df.columns:
        print("Required columns not found")
        return None

    # Extract ontology prefix
    df = metadata_df.copy()
    df['ontology_prefix'] = df['ontology_curie'].apply(
        lambda x: x.split(':')[0] if pd.notna(x) and ':' in str(x) else None
    )

    # Create coverage matrix
    coverage = pd.crosstab(
        df['output_type'],
        df['ontology_prefix'],
        margins=True
    )

    print(coverage)

    return coverage


def main():
    """Main execution function."""
    print("AlphaGenome Tissue Ontology Mapping Tutorial")
    print("=" * 80)

    # Initialize model
    print("\nInitializing AlphaGenome DNA model...")
    dna_model = dna_client.create(API_KEY)

    # 1. Explore output metadata
    metadata_results = explore_output_metadata(dna_model)

    # Save track counts
    track_counts_file = RESULTS_DIR / "track_counts.csv"
    metadata_results['track_counts'].to_csv(track_counts_file)
    print(f"\nSaved track counts to: {track_counts_file}")

    # 2. Extract ontology terms
    human_ontologies = extract_ontology_terms(metadata_results['human_metadata'])
    mouse_ontologies = extract_ontology_terms(metadata_results['mouse_metadata'])

    # 3. Search for specific tissues/cells
    search_terms = ['brain', 'liver', 'heart', 'lung', 'T cell', 'neuron']
    search_results = search_tissue_terms(
        metadata_results['human_metadata'],
        search_terms
    )

    # 4. Analyze ontology coverage
    coverage = analyze_ontology_coverage(metadata_results['human_metadata'])

    if coverage is not None:
        coverage_file = RESULTS_DIR / "ontology_coverage.csv"
        coverage.to_csv(coverage_file)
        print(f"\nSaved ontology coverage to: {coverage_file}")

    # 5. Save all ontology terms to JSON
    ontology_data = {
        'human_ontologies': {k: v for k, v in human_ontologies.items()},
        'mouse_ontologies': {k: v for k, v in mouse_ontologies.items()},
        'search_results': {
            term: {
                'count': result['count'],
                'sample_entries': result['matches'].head(3).to_dict('records') if len(result['matches']) > 0 else []
            }
            for term, result in search_results.items()
        }
    }

    ontology_json_file = RESULTS_DIR / "ontology_terms.json"
    with open(ontology_json_file, 'w') as f:
        json.dump(ontology_data, f, indent=2, default=str)

    print(f"\nSaved ontology terms to: {ontology_json_file}")

    # 6. Create summary of available ontologies
    summary = {
        'total_human_tracks': len(metadata_results['human_metadata']),
        'total_mouse_tracks': len(metadata_results['mouse_metadata']),
        'human_ontology_types': list(human_ontologies.keys()),
        'mouse_ontology_types': list(mouse_ontologies.keys()),
        'human_ontology_counts': {k: len(v) for k, v in human_ontologies.items()},
        'mouse_ontology_counts': {k: len(v) for k, v in mouse_ontologies.items()},
    }

    summary_file = RESULTS_DIR / "ontology_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nSaved summary to: {summary_file}")

    # Print final summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\nTotal human tracks: {summary['total_human_tracks']}")
    print(f"Total mouse tracks: {summary['total_mouse_tracks']}")
    print(f"\nHuman ontology types: {summary['human_ontology_types']}")
    print(f"Mouse ontology types: {summary['mouse_ontology_types']}")
    print(f"\nHuman ontology term counts:")
    for ont, count in summary['human_ontology_counts'].items():
        print(f"  {ont}: {count} terms")

    print("\n" + "=" * 80)
    print("TUTORIAL COMPLETE")
    print("=" * 80)
    print(f"\nResults saved to: {RESULTS_DIR}")
    print("\nGenerated files:")
    print(f"  - track_counts.csv")
    print(f"  - ontology_coverage.csv")
    print(f"  - ontology_terms.json")
    print(f"  - ontology_summary.json")


if __name__ == "__main__":
    main()
