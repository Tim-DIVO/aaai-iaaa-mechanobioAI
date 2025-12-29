#!/usr/bin/env python3
"""
Add impact scores to the domain variants CSV.

Impact score calculation:
1. If both SIFT and PolyPhen available: average of (1-SIFT) and PolyPhen
2. If only SIFT: 1 - SIFT
3. If only PolyPhen: PolyPhen
4. If neither: sample from Beta distribution based on consequence type
"""

import pandas as pd
import numpy as np
from pathlib import Path


def calculate_impact_score(row, random_state=None):
    """
    Calculate impact score from SIFT and PolyPhen scores.
    
    Args:
        row: DataFrame row with sift_score, polyphen_score, consequence_type
        random_state: numpy RandomState for reproducible sampling
    
    Returns:
        tuple: (impact_score, source)
    """
    sift = row['sift_score']
    polyphen = row['polyphen_score']
    consequence = row['consequence_type']
    
    # Check if values are available (not NaN)
    has_sift = pd.notna(sift)
    has_polyphen = pd.notna(polyphen)
    
    if has_sift and has_polyphen:
        # Both available: average
        impact_from_sift = 1 - sift
        impact_from_polyphen = polyphen
        return (impact_from_sift + impact_from_polyphen) / 2, 'both'
    
    elif has_sift:
        # Only SIFT
        return 1 - sift, 'sift_only'
    
    elif has_polyphen:
        # Only PolyPhen
        return polyphen, 'polyphen_only'
    
    else:
        # Neither: sample based on consequence type
        # Beta distribution parameters (a, b) chosen to match expected impact
        consequence_distributions = {
            'stop gained': (8, 2),        # mean ~0.8, high impact
            'frameshift': (8, 2),          # mean ~0.8, high impact
            'inframe deletion': (6, 3),    # mean ~0.67, moderate-high impact
            'delins': (6, 3),              # mean ~0.67, moderate-high impact
            'initiator codon variant': (8, 2),  # mean ~0.8, high impact
            'missense': (3, 3),            # mean ~0.5, variable impact
            'insertion': (3, 3),           # mean ~0.5, variable impact
            '-': (2, 2),                   # mean ~0.5, unknown
        }
        
        # Get distribution parameters, default to (2, 2) for unknown types
        a, b = consequence_distributions.get(consequence, (2, 2))
        
        # Sample from Beta distribution
        if random_state is not None:
            score = random_state.beta(a, b)
        else:
            score = np.random.beta(a, b)
        
        return score, 'sampled'


def add_impact_scores(input_csv: Path, output_csv: Path, seed: int = 42):
    """
    Add impact scores to variant CSV.
    
    Args:
        input_csv: Input CSV file path
        output_csv: Output CSV file path
        seed: Random seed for reproducibility
    """
    print(f"Loading variants from {input_csv}...")
    df = pd.read_csv(input_csv)
    
    print(f"  Total variants: {len(df)}")
    print(f"  Columns: {list(df.columns)}")
    
    # Initialize random state for reproducibility
    random_state = np.random.RandomState(seed)
    
    # Calculate impact scores
    print(f"\nCalculating impact scores (seed={seed})...")
    results = df.apply(lambda row: calculate_impact_score(row, random_state), axis=1)
    df['impact_score'] = results.apply(lambda x: x[0])
    df['impact_source'] = results.apply(lambda x: x[1])
    
    # Summary statistics
    print("\n" + "=" * 80)
    print("IMPACT SCORE SUMMARY")
    print("=" * 80)
    
    source_counts = df['impact_source'].value_counts()
    for source in ['both', 'sift_only', 'polyphen_only', 'sampled']:
        if source in source_counts.index:
            count = source_counts[source]
            pct = 100 * count / len(df)
            print(f"  {source:15s}: {count:5d} ({pct:5.1f}%)")
    
    print(f"\n  Overall statistics:")
    print(f"    Mean:   {df['impact_score'].mean():.3f}")
    print(f"    Median: {df['impact_score'].median():.3f}")
    print(f"    Std:    {df['impact_score'].std():.3f}")
    print(f"    Min:    {df['impact_score'].min():.3f}")
    print(f"    Max:    {df['impact_score'].max():.3f}")
    
    # Statistics by source
    print("\n" + "=" * 80)
    print("STATISTICS BY SOURCE")
    print("=" * 80)
    for source in ['both', 'sift_only', 'polyphen_only', 'sampled']:
        subset = df[df['impact_source'] == source]['impact_score']
        if len(subset) > 0:
            print(f"\n  {source} (n={len(subset)}):")
            print(f"    Mean:   {subset.mean():.3f}")
            print(f"    Median: {subset.median():.3f}")
            print(f"    Std:    {subset.std():.3f}")
    
    # Statistics for sampled by consequence type
    sampled_df = df[df['impact_source'] == 'sampled']
    if len(sampled_df) > 0:
        print("\n" + "=" * 80)
        print("SAMPLED SCORES BY CONSEQUENCE TYPE")
        print("=" * 80)
        for consequence in sampled_df['consequence_type'].unique():
            subset = sampled_df[sampled_df['consequence_type'] == consequence]['impact_score']
            print(f"\n  {consequence} (n={len(subset)}):")
            print(f"    Mean:   {subset.mean():.3f}")
            print(f"    Median: {subset.median():.3f}")
    
    # Write output
    print(f"\n{'=' * 80}")
    print(f"Writing output to {output_csv}...")
    print(f"{'=' * 80}")
    df.to_csv(output_csv, index=False)
    
    print(f"\n✅ Impact scores added successfully!")
    print(f"   Output: {output_csv}")
    print(f"   Total rows: {len(df)}")
    print(f"   New columns: impact_score, impact_source")
    print(f"{'=' * 80}")
    
    return df


if __name__ == '__main__':
    input_file = Path('../protein_csvs/domain_variants.csv')
    output_file = Path('../protein_csvs/domain_variants_with_impact.csv')

    if not input_file.exists():
        print(f"Error: {input_file} not found!")
        exit(1)

    add_impact_scores(input_file, output_file, seed=42)

