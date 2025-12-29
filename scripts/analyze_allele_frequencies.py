#!/usr/bin/env python3
"""
Analyze allele frequency distribution to determine appropriate sampling strategy.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats


def analyze_allele_frequencies(csv_file: Path):
    """Analyze and visualize allele frequency distributions."""
    
    print(f"Loading data from {csv_file}...")
    df = pd.read_csv(csv_file)
    
    print(f"Total variants: {len(df)}")
    
    # Check availability
    has_freq = df['allele_frequency'].notna()
    n_with_freq = has_freq.sum()
    n_without_freq = (~has_freq).sum()
    
    print("\n" + "=" * 80)
    print("ALLELE FREQUENCY AVAILABILITY")
    print("=" * 80)
    print(f"  With frequency:    {n_with_freq:5d} ({100*n_with_freq/len(df):5.1f}%)")
    print(f"  Without frequency: {n_without_freq:5d} ({100*n_without_freq/len(df):5.1f}%)")
    
    # Get variants with frequencies
    df_with_freq = df[has_freq].copy()
    freqs = df_with_freq['allele_frequency'].values

    if len(freqs) == 0:
        print("\n⚠️  No allele frequencies available!")
        return

    # Check for zeros
    n_zeros = (freqs == 0).sum()
    n_nonzero = (freqs > 0).sum()

    print(f"\n  Zero frequencies: {n_zeros} ({100*n_zeros/len(freqs):.1f}%)")
    print(f"  Non-zero frequencies: {n_nonzero} ({100*n_nonzero/len(freqs):.1f}%)")

    # Statistics
    print("\n" + "=" * 80)
    print("ALLELE FREQUENCY STATISTICS (where available)")
    print("=" * 80)
    print(f"  Count:  {len(freqs)}")
    print(f"  Mean:   {freqs.mean():.6f}")
    print(f"  Median: {np.median(freqs):.6f}")
    print(f"  Std:    {freqs.std():.6f}")
    print(f"  Min:    {freqs.min():.6e}")
    print(f"  Max:    {freqs.max():.6f}")
    print(f"\n  Percentiles:")
    for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
        print(f"    {p:2d}%: {np.percentile(freqs, p):.6e}")

    # Filter out zeros for log-transform
    freqs_nonzero = freqs[freqs > 0]

    if len(freqs_nonzero) == 0:
        print("\n⚠️  All frequencies are zero!")
        return

    # Log-transform for better visualization
    log_freqs = np.log10(freqs_nonzero)
    
    print("\n" + "=" * 80)
    print(f"LOG10(ALLELE FREQUENCY) STATISTICS (non-zero only, n={len(freqs_nonzero)})")
    print("=" * 80)
    print(f"  Mean:   {log_freqs.mean():.3f}")
    print(f"  Median: {np.median(log_freqs):.3f}")
    print(f"  Std:    {log_freqs.std():.3f}")
    print(f"  Min:    {log_freqs.min():.3f}")
    print(f"  Max:    {log_freqs.max():.3f}")
    
    # Test for distribution fit
    print("\n" + "=" * 80)
    print("DISTRIBUTION FITTING")
    print("=" * 80)
    
    # Fit normal to log-frequencies
    mu, sigma = stats.norm.fit(log_freqs)
    print(f"\nLog-Normal fit (on log10 scale):")
    print(f"  μ (mean of log10):  {mu:.3f}")
    print(f"  σ (std of log10):   {sigma:.3f}")
    print(f"  → Typical range: 10^{mu-2*sigma:.1f} to 10^{mu+2*sigma:.1f}")
    
    # Fit Beta distribution to original frequencies (non-zero only)
    # Beta needs values in [0, 1], check if we need to scale
    if freqs_nonzero.max() <= 1.0:
        a, b, loc, scale = stats.beta.fit(freqs_nonzero, floc=0, fscale=1)
        print(f"\nBeta distribution fit:")
        print(f"  α (a): {a:.3f}")
        print(f"  β (b): {b:.3f}")
        print(f"  Mean:  {a/(a+b):.6f}")
    else:
        print(f"\nBeta distribution fit: SKIPPED")
        print(f"  (Max frequency {freqs_nonzero.max():.3f} > 1.0, Beta not applicable)")
        a, b = None, None
    
    # Create visualizations
    print("\n" + "=" * 80)
    print("Creating visualizations...")
    print("=" * 80)
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Allele Frequency Distribution Analysis', fontsize=16, fontweight='bold')
    
    # 1. Histogram of raw frequencies (non-zero only for better visualization)
    ax = axes[0, 0]
    ax.hist(freqs_nonzero, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax.axvline(freqs_nonzero.mean(), color='red', linestyle='--', linewidth=2,
               label=f'Mean: {freqs_nonzero.mean():.6f}')
    ax.axvline(np.median(freqs_nonzero), color='orange', linestyle='--', linewidth=2,
               label=f'Median: {np.median(freqs_nonzero):.6f}')
    ax.set_xlabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title(f'Raw Frequencies (non-zero, n={len(freqs_nonzero)})', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 2. Histogram of log10 frequencies
    ax = axes[0, 1]
    ax.hist(log_freqs, bins=50, edgecolor='black', alpha=0.7, color='green')
    ax.axvline(log_freqs.mean(), color='red', linestyle='--', linewidth=2,
               label=f'Mean: {log_freqs.mean():.2f}')
    ax.axvline(np.median(log_freqs), color='orange', linestyle='--', linewidth=2,
               label=f'Median: {np.median(log_freqs):.2f}')
    ax.set_xlabel('Log10(Allele Frequency)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax.set_title('Log10-Transformed Frequencies', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 3. Q-Q plot for log-normal
    ax = axes[0, 2]
    stats.probplot(log_freqs, dist="norm", plot=ax)
    ax.set_title('Q-Q Plot (Log-Normal)', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # 4. Cumulative distribution (non-zero only)
    ax = axes[1, 0]
    sorted_freqs = np.sort(freqs_nonzero)
    cumulative = np.arange(1, len(sorted_freqs) + 1) / len(sorted_freqs)
    ax.plot(sorted_freqs, cumulative, linewidth=2, color='steelblue', label='Empirical')
    ax.set_xlabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax.set_ylabel('Cumulative Probability', fontsize=11, fontweight='bold')
    ax.set_title('Cumulative Distribution (non-zero)', fontsize=12, fontweight='bold')
    ax.set_xscale('log')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 5. Box plot (non-zero only)
    ax = axes[1, 1]
    bp = ax.boxplot([freqs_nonzero], vert=True, patch_artist=True, labels=['Non-zero'])
    bp['boxes'][0].set_facecolor('steelblue')
    bp['boxes'][0].set_alpha(0.6)
    ax.set_ylabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax.set_title('Distribution Summary', fontsize=12, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(alpha=0.3, axis='y')
    
    # 6. Frequency by protein
    ax = axes[1, 2]
    proteins = df_with_freq['protein'].unique()
    freq_by_protein = [df_with_freq[df_with_freq['protein'] == p]['allele_frequency'].values 
                       for p in proteins]
    bp = ax.boxplot(freq_by_protein, tick_labels=proteins, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('steelblue')
        patch.set_alpha(0.6)
    ax.set_ylabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax.set_title('Frequency by Protein', fontsize=12, fontweight='bold')
    ax.set_yscale('log')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()

    output_file = '../plots/allele_frequency_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Plot saved to: {output_file}")
    
    # Recommendations
    print("\n" + "=" * 80)
    print("SAMPLING RECOMMENDATIONS")
    print("=" * 80)
    print(f"\nBased on the analysis, for missing allele frequencies:")
    print(f"\n  Option 1: Log-Normal Distribution")
    print(f"    Sample from: LogNormal(μ={mu:.2f}, σ={sigma:.2f})")
    print(f"    In code: 10 ** np.random.normal({mu:.2f}, {sigma:.2f})")
    print(f"    Typical range: {10**(mu-2*sigma):.2e} to {10**(mu+2*sigma):.2e}")
    
    if a is not None and b is not None:
        print(f"\n  Option 2: Beta Distribution")
        print(f"    Sample from: Beta(α={a:.2f}, β={b:.2f})")
        print(f"    Mean: {a/(a+b):.6f}")
    else:
        print(f"\n  Option 2: Beta Distribution - Not applicable (frequencies > 1)")
    
    print(f"\n  Option 3: Empirical Resampling")
    print(f"    Randomly sample from the {len(freqs)} observed frequencies")
    
    print(f"\n  Recommendation: Use Log-Normal (Option 1)")
    print(f"    - Best fit for rare variant frequencies")
    print(f"    - Captures the heavy tail toward rare variants")
    print(f"    - Biologically realistic")
    print(f"\n  Note: {n_zeros} variants have frequency = 0")
    print(f"    - These are likely singletons or very rare variants")
    print(f"    - For sampling, treat as extremely rare: use minimum observed non-zero value")

    return mu, sigma, a, b, freqs_nonzero.min()


if __name__ == '__main__':
    csv_file = Path('../protein_csvs/domain_variants_with_impact.csv')

    if not csv_file.exists():
        print(f"Error: {csv_file} not found!")
        exit(1)

    analyze_allele_frequencies(csv_file)

