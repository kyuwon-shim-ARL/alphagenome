# AlphaGenome Essential Commands - Execution Results

## Overview
This directory contains the execution results of the AlphaGenome Essential Commands tutorial from `tutorials/essential_commands.ipynb`.

## Files
- `run_essential_commands.py`: Executable script extracting core functionality from the tutorial
- `results.json`: Complete test results in JSON format

## Test Coverage

### 1. Interval Operations (`genome.Interval`)
- ✓ Creating intervals with chromosome, start, end
- ✓ Accessing properties: center(), width
- ✓ Resizing intervals around center point
- ✓ Comparing intervals: overlaps(), contains(), intersect()

**Example Results:**
```
Interval: chr1:1000-1010
Center: 1005
Width: 10
Resized (100bp): chr1:955-1055
Intersection: chr1:1005-1010
```

### 2. Variant Operations (`genome.Variant`)
- ✓ Creating SNVs (single nucleotide variants)
- ✓ Creating insertions
- ✓ Creating deletions (indels)
- ✓ Accessing reference_interval
- ✓ Testing overlap with intervals (reference_overlaps, alternate_overlaps)

**Example Results:**
```
SNV: chr3:10000:A>C
Insertion: chr3:10000:T>CGTCAAT
Deletion: chr3:10000:AGGGATC>C
Reference interval: chr3:9999-10000
```

### 3. TrackData Operations (`track_data.TrackData`)
- ✓ Creating TrackData from numpy arrays and pandas metadata
- ✓ Setting resolution and genomic interval
- ✓ Changing resolution (downsampling and upsampling)
- ✓ Filtering by strand (positive, negative, unstranded)
- ✓ Resizing (cropping and padding)
- ✓ Slicing by position and interval
- ✓ Selecting tracks by name and index
- ✓ Reverse complement transformation

**Example Results:**
```
Shape: [4, 3] (4 positions x 3 tracks)
Positive strand tracks: ['track1']
Negative strand tracks: ['track1']
Unstranded tracks: ['track2']
```

**Resolution Conversion:**
- Downsampling (1bp → 2bp): Sums adjacent values
- Upsampling (2bp → 1bp): Repeats values while preserving sum

## Key Insights

### Indexing Convention
- `genome.Interval`: 0-based indexing (includes start, excludes end)
- `genome.Variant`: 1-based position (compatible with VCF format)

### TrackData Structure
- `values`: numpy.ndarray of shape (sequence_length, num_tracks)
- `metadata`: pandas.DataFrame with track annotations (name, strand, etc.)
- `resolution`: Base pairs per value (1 = single base resolution)
- `interval`: Genomic coordinates of the data

### Strand Awareness
- Tracks can be stranded (+, -) or unstranded (.)
- Reverse complement automatically handles strand switching
- Filtering methods available for each strand type

## Usage
```bash
# Run the essential commands test
uv run python results/essential_commands/run_essential_commands.py

# View results
cat results/essential_commands/results.json
```

## Notes
- All operations are tested without requiring API authentication
- Tests use synthetic data to demonstrate functionality
- Real predictions would require `dna_client.DnaClient` with API key
