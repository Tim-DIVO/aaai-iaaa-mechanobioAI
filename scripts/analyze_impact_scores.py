#!/usr/bin/env python3
"""
Analyze the distribution of impact scores calculated from SIFT and PolyPhen scores.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def calculate_impact_score(row):
    """Calculate impact score from SIFT and PolyPhen scores."""
    sift = row['sift_score']
    polyphen = row['polyphen_score']
    
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
        # Neither
        return np.nan, 'neither'


def analyze_impact_scores(csv_file: Path):
    """Analyze and plot impact score distributions."""
    
    print(f"Loading data from {csv_file}...")
    df = pd.read_csv(csv_file)
    
    print(f"Total variants: {len(df)}")
    
    # Calculate impact scores
    print("\nCalculating impact scores...")
    df[['impact_score', 'score_source']] = df.apply(
        calculate_impact_score, axis=1, result_type='expand'
    )
    
    # Summary statistics
    print("\n" + "=" * 80)
    print("IMPACT SCORE AVAILABILITY")
    print("=" * 80)
    source_counts = df['score_source'].value_counts()
    for source, count in source_counts.items():
        pct = 100 * count / len(df)
        print(f"  {source:15s}: {count:5d} ({pct:5.1f}%)")
    
    # Get variants with impact scores
    df_with_scores = df[df['impact_score'].notna()].copy()
    print(f"\nTotal variants with impact scores: {len(df_with_scores)}")
    
    # Statistics by source
    print("\n" + "=" * 80)
    print("IMPACT SCORE STATISTICS BY SOURCE")
    print("=" * 80)
    for source in ['both', 'sift_only', 'polyphen_only']:
        subset = df_with_scores[df_with_scores['score_source'] == source]['impact_score']
        if len(subset) > 0:
            print(f"\n{source}:")
            print(f"  Count: {len(subset)}")
            print(f"  Mean:  {subset.mean():.3f}")
            print(f"  Std:   {subset.std():.3f}")
            print(f"  Min:   {subset.min():.3f}")
            print(f"  25%:   {subset.quantile(0.25):.3f}")
            print(f"  50%:   {subset.median():.3f}")
            print(f"  75%:   {subset.quantile(0.75):.3f}")
            print(f"  Max:   {subset.max():.3f}")
    
    # Overall statistics
    print("\n" + "=" * 80)
    print("OVERALL IMPACT SCORE STATISTICS")
    print("=" * 80)
    print(f"  Count: {len(df_with_scores)}")
    print(f"  Mean:  {df_with_scores['impact_score'].mean():.3f}")
    print(f"  Std:   {df_with_scores['impact_score'].std():.3f}")
    print(f"  Min:   {df_with_scores['impact_score'].min():.3f}")
    print(f"  25%:   {df_with_scores['impact_score'].quantile(0.25):.3f}")
    print(f"  50%:   {df_with_scores['impact_score'].median():.3f}")
    print(f"  75%:   {df_with_scores['impact_score'].quantile(0.75):.3f}")
    print(f"  Max:   {df_with_scores['impact_score'].max():.3f}")
    
    # Create visualizations
    print("\n" + "=" * 80)
    print("Creating visualizations...")
    print("=" * 80)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Impact Score Distribution Analysis', fontsize=16, fontweight='bold')
    
    # 1. Overall distribution
    ax = axes[0, 0]
    ax.hist(df_with_scores['impact_score'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax.axvline(df_with_scores['impact_score'].mean(), color='red', linestyle='--', 
               linewidth=2, label=f'Mean: {df_with_scores["impact_score"].mean():.3f}')
    ax.axvline(df_with_scores['impact_score'].median(), color='orange', linestyle='--', 
               linewidth=2, label=f'Median: {df_with_scores["impact_score"].median():.3f}')
    ax.set_xlabel('Impact Score', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Overall Impact Score Distribution', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 2. Distribution by source
    ax = axes[0, 1]
    for source, color in [('both', 'green'), ('sift_only', 'blue'), ('polyphen_only', 'orange')]:
        subset = df_with_scores[df_with_scores['score_source'] == source]['impact_score']
        if len(subset) > 0:
            ax.hist(subset, bins=30, alpha=0.5, label=f'{source} (n={len(subset)})', color=color, edgecolor='black')
    ax.set_xlabel('Impact Score', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Impact Score by Source', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    # 3. Box plot by source
    ax = axes[1, 0]
    source_order = ['both', 'sift_only', 'polyphen_only']
    data_for_box = [df_with_scores[df_with_scores['score_source'] == s]['impact_score'].values 
                    for s in source_order]
    bp = ax.boxplot(data_for_box, labels=source_order, patch_artist=True)
    for patch, color in zip(bp['boxes'], ['green', 'blue', 'orange']):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_ylabel('Impact Score', fontsize=12)
    ax.set_title('Impact Score Distribution by Source', fontsize=13, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    # 4. Cumulative distribution
    ax = axes[1, 1]
    sorted_scores = np.sort(df_with_scores['impact_score'])
    cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    ax.plot(sorted_scores, cumulative, linewidth=2, color='steelblue')
    ax.axhline(0.5, color='red', linestyle='--', alpha=0.5, label='50th percentile')
    ax.axhline(0.25, color='orange', linestyle='--', alpha=0.5, label='25th percentile')
    ax.axhline(0.75, color='orange', linestyle='--', alpha=0.5, label='75th percentile')
    ax.set_xlabel('Impact Score', fontsize=12)
    ax.set_ylabel('Cumulative Probability', fontsize=12)
    ax.set_title('Cumulative Distribution Function', fontsize=13, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()

    output_file = '../plots/impact_score_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Plot saved to: {output_file}")

    # Save the dataframe with impact scores for inspection
    output_csv = '../protein_csvs/domain_variants_with_impact.csv'
    df.to_csv(output_csv, index=False)
    print(f"✅ Data with impact scores saved to: {output_csv}")
    
    return df


if __name__ == '__main__':
    csv_file = Path('../protein_csvs/domain_variants.csv')

    if not csv_file.exists():
        print(f"Error: {csv_file} not found!")
        exit(1)

    analyze_impact_scores(csv_file)

