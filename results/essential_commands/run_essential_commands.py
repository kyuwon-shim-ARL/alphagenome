"""
Essential Commands Tutorial - Code Execution
Extracted from tutorials/essential_commands.ipynb
"""

import json
import numpy as np
import pandas as pd
from alphagenome.data import genome, track_data

def test_interval_operations():
    """Test genome.Interval operations"""
    results = {}

    # Create interval
    interval = genome.Interval(chromosome='chr1', start=1_000, end=1_010)
    results['interval_created'] = str(interval)

    # Interval properties
    results['interval_center'] = interval.center()
    results['interval_width'] = interval.width

    # Resize
    resized = interval.resize(100)
    results['interval_resized'] = str(resized)

    # Compare intervals
    second_interval = genome.Interval(chromosome='chr1', start=1_005, end=1_015)
    results['overlaps'] = interval.overlaps(second_interval)
    results['contains'] = interval.contains(second_interval)
    results['intersect'] = str(interval.intersect(second_interval))

    return results

def test_variant_operations():
    """Test genome.Variant operations"""
    results = {}

    # Create variant
    variant = genome.Variant(
        chromosome='chr3', position=10_000, reference_bases='A', alternate_bases='C'
    )
    results['variant_created'] = str(variant)

    # Insertion variant
    insertion = genome.Variant(
        chromosome='chr3',
        position=10_000,
        reference_bases='T',
        alternate_bases='CGTCAAT',
    )
    results['insertion_created'] = str(insertion)

    # Deletion variant
    deletion = genome.Variant(
        chromosome='chr3',
        position=10_000,
        reference_bases='AGGGATC',
        alternate_bases='C',
    )
    results['deletion_created'] = str(deletion)

    # Reference interval
    variant = genome.Variant(
        chromosome='chr3', position=10_000, reference_bases='A', alternate_bases='T'
    )
    results['reference_interval'] = str(variant.reference_interval)

    # Overlap with interval
    variant = genome.Variant(
        chromosome='chr3',
        position=10_000,
        reference_bases='T',
        alternate_bases='CGTCAAT',
    )

    interval = genome.Interval(chromosome='chr3', start=10_005, end=10_010)
    results['reference_overlaps'] = variant.reference_overlaps(interval)
    results['alternate_overlaps'] = variant.alternate_overlaps(interval)

    return results

def test_trackdata_operations():
    """Test TrackData operations"""
    results = {}

    # Create TrackData from scratch
    values = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]).astype(np.float32)
    metadata = pd.DataFrame({
        'name': ['track1', 'track1', 'track2'],
        'strand': ['+', '-', '.'],
    })

    tdata = track_data.TrackData(values=values, metadata=metadata)
    results['trackdata_created'] = True
    results['trackdata_shape'] = list(tdata.values.shape)

    # With resolution
    interval = genome.Interval(chromosome='chr1', start=1_000, end=1_004)
    tdata = track_data.TrackData(
        values=values, metadata=metadata, resolution=1, interval=interval
    )
    results['trackdata_with_resolution'] = True

    # Change resolution - downsample
    tdata_downsampled = tdata.change_resolution(resolution=2)
    results['downsampled_values'] = tdata_downsampled.values.tolist()

    # Change resolution - upsample
    tdata_upsampled = tdata_downsampled.change_resolution(resolution=1)
    results['upsampled_values'] = tdata_upsampled.values.tolist()

    # Filtering
    results['positive_strand_tracks'] = tdata.filter_to_positive_strand().metadata.name.values.tolist()
    results['negative_strand_tracks'] = tdata.filter_to_negative_strand().metadata.name.values.tolist()
    results['unstranded_tracks'] = tdata.filter_to_unstranded().metadata.name.values.tolist()

    # Resizing - smaller (cropping)
    resized_small = tdata.resize(width=2)
    results['resized_small'] = resized_small.values.tolist()

    # Resizing - bigger (padding)
    resized_big = tdata.resize(width=8)
    results['resized_big'] = resized_big.values.tolist()

    # Slicing by position
    sliced_pos = tdata.slice_by_positions(start=2, end=4)
    results['sliced_by_position'] = sliced_pos.values.tolist()

    # Slicing by interval
    sliced_int = tdata.slice_by_interval(
        genome.Interval(chromosome='chr1', start=1_002, end=1_004)
    )
    results['sliced_by_interval'] = sliced_int.values.tolist()

    # Subset tracks by name
    track1_tdata = tdata.select_tracks_by_name(names='track1')
    results['track1_values'] = track1_tdata.values.tolist()
    results['track1_metadata'] = track1_tdata.metadata.name.values.tolist()

    # Select tracks by index
    indexed_tdata = tdata.select_tracks_by_index(idx=[0, 1])
    results['tracks_by_index'] = indexed_tdata.values.tolist()

    # Reverse complement
    interval_stranded = genome.Interval(
        chromosome='chr1', start=1_000, end=1_004, strand='+'
    )
    tdata_stranded = track_data.TrackData(
        values=values, metadata=metadata, resolution=1, interval=interval_stranded
    )
    results['reverse_complement'] = tdata_stranded.reverse_complement().values.tolist()

    return results

def main():
    """Run all tests and save results"""
    all_results = {}

    print("Testing Interval operations...")
    all_results['interval_operations'] = test_interval_operations()
    print("✓ Interval operations completed")

    print("\nTesting Variant operations...")
    all_results['variant_operations'] = test_variant_operations()
    print("✓ Variant operations completed")

    print("\nTesting TrackData operations...")
    all_results['trackdata_operations'] = test_trackdata_operations()
    print("✓ TrackData operations completed")

    # Save results
    output_path = '/home/kyuwon/projects/alphagenome/results/essential_commands/results.json'
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f"\n✓ Results saved to {output_path}")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY OF ESSENTIAL COMMANDS")
    print("="*60)

    print("\n1. Interval Operations:")
    print(f"   - Created: {all_results['interval_operations']['interval_created']}")
    print(f"   - Center: {all_results['interval_operations']['interval_center']}")
    print(f"   - Width: {all_results['interval_operations']['interval_width']}")
    print(f"   - Overlaps: {all_results['interval_operations']['overlaps']}")

    print("\n2. Variant Operations:")
    print(f"   - Variant: {all_results['variant_operations']['variant_created']}")
    print(f"   - Reference overlaps: {all_results['variant_operations']['reference_overlaps']}")
    print(f"   - Alternate overlaps: {all_results['variant_operations']['alternate_overlaps']}")

    print("\n3. TrackData Operations:")
    print(f"   - Shape: {all_results['trackdata_operations']['trackdata_shape']}")
    print(f"   - Positive strand tracks: {all_results['trackdata_operations']['positive_strand_tracks']}")
    print(f"   - Negative strand tracks: {all_results['trackdata_operations']['negative_strand_tracks']}")
    print(f"   - Unstranded tracks: {all_results['trackdata_operations']['unstranded_tracks']}")

    print("\n" + "="*60)

if __name__ == '__main__':
    main()
