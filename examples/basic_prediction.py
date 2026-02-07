#!/usr/bin/env python3
"""
Basic AlphaGenome usage example.

This script demonstrates how to:
1. Connect to the AlphaGenome API
2. Define a genomic interval
3. Define a variant
4. Run a prediction
5. Visualize results

Before running, ensure ALPHAGENOME_API_KEY is set in the .env file.

Note: Import paths are based on available documentation. If imports fail,
verify against official docs at https://www.alphagenomedocs.com/
"""

import os
import sys

# Load environment variables from .env file
from dotenv import load_dotenv

# Load .env from project root
env_path = os.path.join(os.path.dirname(__file__), '..', '.env')
load_dotenv(env_path)

def get_api_key():
    """Get API key from environment."""
    api_key = os.environ.get('ALPHAGENOME_API_KEY')

    if not api_key or api_key == '<your_api_key_here>' or api_key == 'your_api_key_here':
        print("ERROR: ALPHAGENOME_API_KEY not configured")
        print("Please set it in your .env file:")
        print("  1. Open /home/kyuwon/projects/alphagenome/.env")
        print("  2. Replace <your_api_key_here> with your actual API key")
        print("Get your key at: https://deepmind.google.com/science/alphagenome")
        sys.exit(1)

    return api_key

def main():
    # Import AlphaGenome modules
    from alphagenome.data import genome
    from alphagenome.models import dna_client
    from alphagenome.visualization import plot_components
    import matplotlib.pyplot as plt

    # Load API key
    api_key = get_api_key()

    # Create model client
    print("Connecting to AlphaGenome API...")
    model = dna_client.create(api_key)
    print("Connected!")

    # Define genomic interval (region of interest on chromosome 22)
    interval = genome.Interval(
        chromosome='chr22',
        start=35677410,
        end=36725986
    )
    print(f"Analyzing interval: {interval.chromosome}:{interval.start}-{interval.end}")

    # Define a variant (SNP: A -> C)
    variant = genome.Variant(
        chromosome='chr22',
        position=36201698,
        reference_bases='A',
        alternate_bases='C',
    )
    print(f"Variant: {variant.chromosome}:{variant.position} {variant.reference_bases}>{variant.alternate_bases}")

    # Run prediction
    print("Running prediction (this may take a moment)...")
    outputs = model.predict_variant(
        interval=interval,
        variant=variant,
        ontology_terms=['UBERON:0001157'],  # Small intestine
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
    )

    print("Prediction complete!")
    print(f"Output type: {type(outputs).__name__}")
    print(f"Output: {outputs}")

    # Optional: Visualize results
    # plot_components.plot_predictions(outputs)
    # plt.savefig('output/prediction_results.png')
    # print("Results saved to output/prediction_results.png")

    return outputs

if __name__ == "__main__":
    main()
