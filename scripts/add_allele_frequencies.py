#!/usr/bin/env python3
"""
Add allele frequencies to variants that are missing them.

Uses Log-Normal distribution based on observed frequency distribution:
- μ = -5.217 (mean of log10 frequencies)
- σ = 0.924 (std of log10 frequencies)
"""

import pandas as pd
import numpy as np
from pathlib import Path


def add_allele_frequencies(input_csv: Path, output_csv: Path, seed: int = 42):
    """
    Add allele frequencies to variants missing them.
    
    Args:
        input_csv: Input CSV file path
        output_csv: Output CSV file path
        seed: Random seed for reproducibility
    """
    print(f"Loading data from {input_csv}...")
    df = pd.read_csv(input_csv)
    
    print(f"  Total variants: {len(df)}")
    
    # Check current status
    has_freq = df['allele_frequency'].notna()
    n_with_freq = has_freq.sum()
    n_without_freq = (~has_freq).sum()
    
    print("\n" + "=" * 80)
    print("CURRENT STATUS")
    print("=" * 80)
    print(f"  With frequency:    {n_with_freq:5d} ({100*n_with_freq/len(df):5.1f}%)")
    print(f"  Without frequency: {n_without_freq:5d} ({100*n_without_freq/len(df):5.1f}%)")
    
    if n_without_freq == 0:
        print("\n✅ All variants already have allele frequencies!")
        return df
    
    # Get observed frequencies for statistics
    observed_freqs = df[has_freq]['allele_frequency'].values
    observed_nonzero = observed_freqs[observed_freqs > 0]
    
    print(f"\n  Observed frequency statistics:")
    print(f"    Mean:   {observed_freqs.mean():.6e}")
    print(f"    Median: {np.median(observed_freqs):.6e}")
    print(f"    Min (non-zero): {observed_nonzero.min():.6e}")
    print(f"    Max:    {observed_freqs.max():.6e}")
    
    # Log-Normal parameters from analysis
    mu = -5.217  # mean of log10(frequency)
    sigma = 0.924  # std of log10(frequency)
    
    print("\n" + "=" * 80)
    print("SAMPLING STRATEGY")
    print("=" * 80)
    print(f"  Distribution: Log-Normal")
    print(f"  Parameters: μ={mu:.3f}, σ={sigma:.3f} (on log10 scale)")
    print(f"  Formula: 10^(N(μ={mu:.3f}, σ={sigma:.3f}))")
    print(f"  Expected range: ~10^{mu-2*sigma:.1f} to ~10^{mu+2*sigma:.1f}")
    
    # Initialize random state
    random_state = np.random.RandomState(seed)
    
    # Sample frequencies for missing values
    print(f"\n  Sampling {n_without_freq} missing frequencies (seed={seed})...")
    
    # Sample from log-normal distribution
    log_freqs_sampled = random_state.normal(mu, sigma, size=n_without_freq)
    freqs_sampled = 10 ** log_freqs_sampled
    
    # Clip to reasonable range (min observed to 1.0)
    min_freq = observed_nonzero.min()
    freqs_sampled = np.clip(freqs_sampled, min_freq, 1.0)
    
    # Add to dataframe
    df.loc[~has_freq, 'allele_frequency'] = freqs_sampled
    
    # Add source column
    df['frequency_source'] = 'observed'
    df.loc[~has_freq, 'frequency_source'] = 'sampled'
    
    # Statistics on sampled frequencies
    print("\n" + "=" * 80)
    print("SAMPLED FREQUENCY STATISTICS")
    print("=" * 80)
    print(f"  Count:  {len(freqs_sampled)}")
    print(f"  Mean:   {freqs_sampled.mean():.6e}")
    print(f"  Median: {np.median(freqs_sampled):.6e}")
    print(f"  Std:    {freqs_sampled.std():.6e}")
    print(f"  Min:    {freqs_sampled.min():.6e}")
    print(f"  Max:    {freqs_sampled.max():.6e}")
    
    # Overall statistics
    all_freqs = df['allele_frequency'].values
    print("\n" + "=" * 80)
    print("OVERALL FREQUENCY STATISTICS (after sampling)")
    print("=" * 80)
    print(f"  Count:  {len(all_freqs)}")
    print(f"  Mean:   {all_freqs.mean():.6e}")
    print(f"  Median: {np.median(all_freqs):.6e}")
    print(f"  Std:    {all_freqs.std():.6e}")
    print(f"  Min:    {all_freqs.min():.6e}")
    print(f"  Max:    {all_freqs.max():.6e}")
    
    # Comparison
    print("\n" + "=" * 80)
    print("COMPARISON: OBSERVED vs SAMPLED")
    print("=" * 80)
    print(f"  {'Metric':<15s} {'Observed':>15s} {'Sampled':>15s} {'Ratio':>10s}")
    print(f"  {'-'*15} {'-'*15} {'-'*15} {'-'*10}")
    print(f"  {'Mean':<15s} {observed_freqs.mean():>15.6e} {freqs_sampled.mean():>15.6e} "
          f"{freqs_sampled.mean()/observed_freqs.mean():>10.2f}")
    print(f"  {'Median':<15s} {np.median(observed_freqs):>15.6e} {np.median(freqs_sampled):>15.6e} "
          f"{np.median(freqs_sampled)/np.median(observed_freqs):>10.2f}")
    
    # Write output
    print(f"\n{'=' * 80}")
    print(f"Writing output to {output_csv}...")
    print(f"{'=' * 80}")
    df.to_csv(output_csv, index=False)
    
    print(f"\n✅ Allele frequencies added successfully!")
    print(f"   Output: {output_csv}")
    print(f"   Total rows: {len(df)}")
    print(f"   Missing frequencies: {df['allele_frequency'].isna().sum()}")
    print(f"   New column: frequency_source (observed/sampled)")
    print(f"{'=' * 80}")
    
    return df


if __name__ == '__main__':
    input_file = Path('../protein_csvs/domain_variants_with_impact.csv')
    output_file = Path('../protein_csvs/domain_variants_complete.csv')

    if not input_file.exists():
        print(f"Error: {input_file} not found!")
        exit(1)

    add_allele_frequencies(input_file, output_file, seed=42)

