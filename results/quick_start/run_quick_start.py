#!/usr/bin/env python3
"""Execute AlphaGenome Quick Start Tutorial"""

import os
import sys
import json
from datetime import datetime

# Setup
os.chdir('/home/kyuwon/projects/alphagenome')
from dotenv import load_dotenv
load_dotenv('.env')

# Imports
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.interpretation import ism
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

results = {"tutorial": "quick_start", "timestamp": datetime.now().isoformat(), "sections": []}

def log_result(section, data):
    results["sections"].append({"section": section, "data": data})
    print(f"[{section}] {data}")

# 1. Create model
print("="*60)
print("AlphaGenome Quick Start Tutorial")
print("="*60)

api_key = os.environ.get('ALPHAGENOME_API_KEY')
dna_model = dna_client.create(api_key)
log_result("Model Created", "Success")

# 2. List output types
output_types = [output.name for output in dna_client.OutputType]
log_result("Output Types", output_types)

# 3. Predict sequence (DNASE for Lung)
print("\n--- Predict Sequence: DNase for Lung ---")
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    requested_outputs=[dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048'],  # Lung
)
log_result("DNase Prediction Shape", str(output.dnase.values.shape))
log_result("DNase Metadata", output.dnase.metadata.to_dict())

# 4. Multiple outputs (CAGE + DNase for Lung + Brain)
print("\n--- Predict Multiple Outputs ---")
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    requested_outputs=[dna_client.OutputType.CAGE, dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048', 'UBERON:0000955'],  # Lung, Brain
)
log_result("DNASE Shape", str(output.dnase.values.shape))
log_result("CAGE Shape", str(output.cage.values.shape))
log_result("CAGE Metadata", output.cage.metadata.to_dict())

# 5. Predict interval (CYP2B6 gene, RNA-seq in liver)
print("\n--- Predict Interval: CYP2B6 RNA-seq ---")
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)
gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
gtf_transcripts = gene_annotation.filter_to_longest_transcript(gtf_transcripts)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)

interval = gene_annotation.get_gene_interval(gtf, gene_symbol='CYP2B6')
log_result("CYP2B6 Interval", str(interval))
interval = interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

output = dna_model.predict_interval(
    interval=interval,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001114'],  # Right liver lobe
)
log_result("RNA-seq Shape", str(output.rna_seq.values.shape))

# Save visualization
transcripts = transcript_extractor.extract(interval)
log_result("Transcripts in Interval", len(transcripts))

plot_components.plot(
    components=[
        plot_components.TranscriptAnnotation(transcripts),
        plot_components.Tracks(output.rna_seq),
    ],
    interval=output.rna_seq.interval,
)
plt.savefig('results/quick_start/cyp2b6_rna_seq.png', dpi=150, bbox_inches='tight')
plt.close()
log_result("Saved Plot", "cyp2b6_rna_seq.png")

# 6. Variant prediction
print("\n--- Variant Prediction ---")
variant = genome.Variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',
    alternate_bases='C',
)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

variant_output = dna_model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0001157'],  # Colon
)
log_result("Variant REF Shape", str(variant_output.reference.rna_seq.values.shape))
log_result("Variant ALT Shape", str(variant_output.alternate.rna_seq.values.shape))

# Save variant plot
transcripts = transcript_extractor.extract(interval)
plot_components.plot(
    [
        plot_components.TranscriptAnnotation(transcripts),
        plot_components.OverlaidTracks(
            tdata={'REF': variant_output.reference.rna_seq, 'ALT': variant_output.alternate.rna_seq},
            colors={'REF': 'dimgrey', 'ALT': 'red'},
        ),
    ],
    interval=variant_output.reference.rna_seq.interval.resize(2**15),
    annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
)
plt.savefig('results/quick_start/variant_effect.png', dpi=150, bbox_inches='tight')
plt.close()
log_result("Saved Plot", "variant_effect.png")

# 7. Variant scoring
print("\n--- Variant Scoring ---")
variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']
variant_scores = dna_model.score_variant(
    interval=interval, variant=variant, variant_scorers=[variant_scorer]
)
log_result("Variant Scores Shape", str(variant_scores[0].X.shape))
log_result("Genes Scored", variant_scores[0].obs['gene_name'].tolist()[:10])

# Tidy scores
tidy = variant_scorers.tidy_scores([variant_scores[0]], match_gene_strand=True)
log_result("Tidy Scores Shape", str(tidy.shape))
tidy.to_csv('results/quick_start/variant_scores.csv', index=False)
log_result("Saved CSV", "variant_scores.csv")

# 8. ISM (smaller example)
print("\n--- In Silico Mutagenesis ---")
sequence_interval = genome.Interval('chr20', 3_753_000, 3_753_400)
sequence_interval = sequence_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)
ism_interval = sequence_interval.resize(64)  # Smaller for speed

dnase_variant_scorer = variant_scorers.CenterMaskScorer(
    requested_output=dna_client.OutputType.DNASE,
    width=501,
    aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
)

ism_scores = dna_model.score_ism_variants(
    interval=sequence_interval,
    ism_interval=ism_interval,
    variant_scorers=[dnase_variant_scorer],
)
log_result("ISM Variants Scored", len(ism_scores))

# 9. Mouse prediction
print("\n--- Mouse Prediction ---")
output = dna_model.predict_sequence(
    sequence='GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N'),
    organism=dna_client.Organism.MUS_MUSCULUS,
    requested_outputs=[dna_client.OutputType.DNASE],
    ontology_terms=['UBERON:0002048'],  # Lung
)
log_result("Mouse DNase Shape", str(output.dnase.values.shape))

# Save results
print("\n" + "="*60)
print("Quick Start Tutorial Complete!")
print("="*60)

with open('results/quick_start/results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("Results saved to results/quick_start/results.json")
