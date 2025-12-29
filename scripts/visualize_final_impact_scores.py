#!/usr/bin/env python3
"""
Create comprehensive visualizations of the final impact scores.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def visualize_impact_scores(csv_file: Path):
    """Create comprehensive visualizations of impact scores."""
    
    print(f"Loading data from {csv_file}...")
    df = pd.read_csv(csv_file)
    
    print(f"Total variants: {len(df)}")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    fig.suptitle('Impact Score Analysis - Complete Dataset', 
                 fontsize=18, fontweight='bold', y=0.995)
    
    # 1. Overall distribution (large plot)
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.hist(df['impact_score'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.axvline(df['impact_score'].mean(), color='red', linestyle='--', 
                linewidth=2, label=f'Mean: {df["impact_score"].mean():.3f}')
    ax1.axvline(df['impact_score'].median(), color='orange', linestyle='--', 
                linewidth=2, label=f'Median: {df["impact_score"].median():.3f}')
    ax1.set_xlabel('Impact Score', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax1.set_title(f'Overall Impact Score Distribution (n={len(df)})', 
                  fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(alpha=0.3)
    
    # 2. Summary statistics table
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.axis('off')
    
    stats_data = []
    for source in ['both', 'sift_only', 'polyphen_only', 'sampled', 'ALL']:
        if source == 'ALL':
            subset = df['impact_score']
            n = len(df)
        else:
            subset = df[df['impact_source'] == source]['impact_score']
            n = len(subset)
        
        if n > 0:
            stats_data.append([
                source,
                f'{n}',
                f'{subset.mean():.3f}',
                f'{subset.median():.3f}',
                f'{subset.std():.3f}'
            ])
    
    table = ax2.table(cellText=stats_data,
                      colLabels=['Source', 'N', 'Mean', 'Median', 'Std'],
                      cellLoc='center',
                      loc='center',
                      bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # Style header
    for i in range(5):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Style ALL row
    for i in range(5):
        table[(len(stats_data), i)].set_facecolor('#E7E6E6')
        table[(len(stats_data), i)].set_text_props(weight='bold')
    
    ax2.set_title('Summary Statistics', fontsize=13, fontweight='bold', pad=20)
    
    # 3. Distribution by source (stacked histogram)
    ax3 = fig.add_subplot(gs[1, :2])
    sources = ['both', 'sift_only', 'polyphen_only', 'sampled']
    colors = ['#2E7D32', '#1976D2', '#F57C00', '#7B1FA2']
    
    data_by_source = [df[df['impact_source'] == s]['impact_score'].values for s in sources]
    ax3.hist(data_by_source, bins=30, label=sources, color=colors, 
             alpha=0.7, edgecolor='black', linewidth=0.5)
    ax3.set_xlabel('Impact Score', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax3.set_title('Impact Score Distribution by Source', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(alpha=0.3)
    
    # 4. Box plot by source
    ax4 = fig.add_subplot(gs[1, 2])
    source_order = ['both', 'sift_only', 'polyphen_only', 'sampled']
    data_for_box = [df[df['impact_source'] == s]['impact_score'].values for s in source_order]
    bp = ax4.boxplot(data_for_box, tick_labels=source_order, patch_artist=True)
    
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    
    ax4.set_ylabel('Impact Score', fontsize=12, fontweight='bold')
    ax4.set_title('Distribution by Source', fontsize=13, fontweight='bold')
    ax4.tick_params(axis='x', rotation=45)
    ax4.grid(alpha=0.3, axis='y')
    
    # 5. Distribution by protein
    ax5 = fig.add_subplot(gs[2, 0])
    protein_means = df.groupby('protein')['impact_score'].mean().sort_values(ascending=False)
    bars = ax5.bar(range(len(protein_means)), protein_means.values, 
                   color='steelblue', edgecolor='black', alpha=0.7)
    ax5.set_xticks(range(len(protein_means)))
    ax5.set_xticklabels(protein_means.index, rotation=45, ha='right')
    ax5.set_ylabel('Mean Impact Score', fontsize=11, fontweight='bold')
    ax5.set_title('Mean Impact Score by Protein', fontsize=12, fontweight='bold')
    ax5.grid(alpha=0.3, axis='y')
    ax5.set_ylim([0, 1])
    
    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, protein_means.values)):
        ax5.text(bar.get_x() + bar.get_width()/2, val + 0.02, f'{val:.2f}',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # 6. Sampled scores by consequence type
    ax6 = fig.add_subplot(gs[2, 1])
    sampled_df = df[df['impact_source'] == 'sampled']
    
    if len(sampled_df) > 0:
        consequence_means = sampled_df.groupby('consequence_type')['impact_score'].mean().sort_values(ascending=False)
        consequence_counts = sampled_df.groupby('consequence_type').size()
        
        bars = ax6.barh(range(len(consequence_means)), consequence_means.values, 
                       color='purple', edgecolor='black', alpha=0.6)
        ax6.set_yticks(range(len(consequence_means)))
        ax6.set_yticklabels([f'{c} (n={consequence_counts[c]})' for c in consequence_means.index], 
                           fontsize=9)
        ax6.set_xlabel('Mean Impact Score', fontsize=11, fontweight='bold')
        ax6.set_title('Sampled Scores by Consequence', fontsize=12, fontweight='bold')
        ax6.grid(alpha=0.3, axis='x')
        ax6.set_xlim([0, 1])
        
        # Add value labels
        for i, (bar, val) in enumerate(zip(bars, consequence_means.values)):
            ax6.text(val + 0.02, bar.get_y() + bar.get_height()/2, f'{val:.2f}',
                    ha='left', va='center', fontsize=8, fontweight='bold')
    
    # 7. Cumulative distribution
    ax7 = fig.add_subplot(gs[2, 2])
    for source, color in zip(sources, colors):
        subset = df[df['impact_source'] == source]['impact_score'].values
        if len(subset) > 0:
            sorted_scores = np.sort(subset)
            cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
            ax7.plot(sorted_scores, cumulative, linewidth=2, color=color, 
                    label=f'{source} (n={len(subset)})', alpha=0.8)
    
    ax7.set_xlabel('Impact Score', fontsize=11, fontweight='bold')
    ax7.set_ylabel('Cumulative Probability', fontsize=11, fontweight='bold')
    ax7.set_title('Cumulative Distribution', fontsize=12, fontweight='bold')
    ax7.legend(fontsize=8)
    ax7.grid(alpha=0.3)
    
    # Save figure
    output_file = '../plots/final_impact_score_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Visualization saved to: {output_file}")
    
    plt.close()


if __name__ == '__main__':
    csv_file = Path('../protein_csvs/domain_variants_with_impact.csv')

    if not csv_file.exists():
        print(f"Error: {csv_file} not found!")
        exit(1)

    visualize_impact_scores(csv_file)
    print("\n✅ All done!")

