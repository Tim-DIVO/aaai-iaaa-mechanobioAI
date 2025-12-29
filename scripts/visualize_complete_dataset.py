#!/usr/bin/env python3
"""
Create comprehensive visualization of the complete dataset with impact scores and allele frequencies.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def visualize_complete_dataset(csv_file: Path):
    """Create comprehensive visualizations of the complete dataset."""
    
    print(f"Loading data from {csv_file}...")
    df = pd.read_csv(csv_file)
    
    print(f"Total variants: {len(df)}")
    print(f"Columns: {list(df.columns)}")
    
    # Create figure
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)
    
    fig.suptitle('Complete Dataset: Impact Scores & Allele Frequencies', 
                 fontsize=18, fontweight='bold', y=0.995)
    
    # 1. Allele frequency distribution (observed vs sampled)
    ax1 = fig.add_subplot(gs[0, 0])
    observed_freqs = df[df['frequency_source'] == 'observed']['allele_frequency'].values
    sampled_freqs = df[df['frequency_source'] == 'sampled']['allele_frequency'].values
    
    # Remove zeros for log plot
    observed_nonzero = observed_freqs[observed_freqs > 0]
    sampled_nonzero = sampled_freqs[sampled_freqs > 0]
    
    ax1.hist([observed_nonzero, sampled_nonzero], bins=40, label=['Observed', 'Sampled'],
             color=['steelblue', 'orange'], alpha=0.6, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax1.set_title('Allele Frequency: Observed vs Sampled', fontsize=12, fontweight='bold')
    ax1.set_xscale('log')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # 2. Impact score distribution (by source)
    ax2 = fig.add_subplot(gs[0, 1])
    for source, color in [('both', '#2E7D32'), ('sift_only', '#1976D2'), 
                          ('polyphen_only', '#F57C00'), ('sampled', '#7B1FA2')]:
        subset = df[df['impact_source'] == source]['impact_score'].values
        if len(subset) > 0:
            ax2.hist(subset, bins=30, alpha=0.5, label=f'{source} (n={len(subset)})',
                    color=color, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Impact Score', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax2.set_title('Impact Score Distribution', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3)
    
    # 3. 2D scatter: Impact vs Frequency
    ax3 = fig.add_subplot(gs[0, 2])
    
    # Color by frequency source
    for freq_src, color, marker in [('observed', 'steelblue', 'o'), ('sampled', 'orange', 's')]:
        subset = df[df['frequency_source'] == freq_src]
        freqs_plot = subset['allele_frequency'].values
        freqs_plot = np.maximum(freqs_plot, 1e-7)  # Avoid log(0)
        ax3.scatter(freqs_plot, subset['impact_score'].values, 
                   alpha=0.4, s=20, color=color, marker=marker, label=freq_src,
                   edgecolors='black', linewidths=0.3)
    
    ax3.set_xlabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Impact Score', fontsize=11, fontweight='bold')
    ax3.set_title('Impact Score vs Allele Frequency', fontsize=12, fontweight='bold')
    ax3.set_xscale('log')
    ax3.legend()
    ax3.grid(alpha=0.3)
    
    # 4. Box plot: Frequency by protein
    ax4 = fig.add_subplot(gs[1, 0])
    proteins = sorted(df['protein'].unique())
    freq_by_protein = [df[df['protein'] == p]['allele_frequency'].values for p in proteins]
    freq_by_protein_nonzero = [[f for f in freqs if f > 0] for freqs in freq_by_protein]
    
    bp = ax4.boxplot(freq_by_protein_nonzero, tick_labels=proteins, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('steelblue')
        patch.set_alpha(0.6)
    ax4.set_ylabel('Allele Frequency', fontsize=11, fontweight='bold')
    ax4.set_title('Frequency Distribution by Protein', fontsize=12, fontweight='bold')
    ax4.set_yscale('log')
    ax4.tick_params(axis='x', rotation=45)
    ax4.grid(alpha=0.3, axis='y')
    
    # 5. Box plot: Impact by protein
    ax5 = fig.add_subplot(gs[1, 1])
    impact_by_protein = [df[df['protein'] == p]['impact_score'].values for p in proteins]
    
    bp = ax5.boxplot(impact_by_protein, tick_labels=proteins, patch_artist=True)
    colors_protein = ['#E57373', '#81C784', '#64B5F6', '#FFD54F', '#BA68C8']
    for patch, color in zip(bp['boxes'], colors_protein):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax5.set_ylabel('Impact Score', fontsize=11, fontweight='bold')
    ax5.set_title('Impact Score by Protein', fontsize=12, fontweight='bold')
    ax5.tick_params(axis='x', rotation=45)
    ax5.grid(alpha=0.3, axis='y')
    
    # 6. Summary statistics table
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')
    
    stats_data = [
        ['Total Variants', f'{len(df)}'],
        ['', ''],
        ['Impact Scores:', ''],
        ['  Mean', f'{df["impact_score"].mean():.3f}'],
        ['  Median', f'{df["impact_score"].median():.3f}'],
        ['', ''],
        ['Allele Freq:', ''],
        ['  Mean', f'{df["allele_frequency"].mean():.2e}'],
        ['  Median', f'{df["allele_frequency"].median():.2e}'],
        ['', ''],
        ['Sources:', ''],
        ['  Observed freq', f'{(df["frequency_source"]=="observed").sum()}'],
        ['  Sampled freq', f'{(df["frequency_source"]=="sampled").sum()}'],
    ]
    
    table = ax6.table(cellText=stats_data, cellLoc='left', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    for i in range(len(stats_data)):
        if stats_data[i][0] in ['Impact Scores:', 'Allele Freq:', 'Sources:']:
            table[(i, 0)].set_facecolor('#E7E6E6')
            table[(i, 0)].set_text_props(weight='bold')
    
    ax6.set_title('Dataset Summary', fontsize=12, fontweight='bold', pad=20)
    
    # 7. Heatmap: Mean impact by protein and domain (top domains)
    ax7 = fig.add_subplot(gs[2, :2])
    
    # Get top domains by variant count
    top_domains = df['domain'].value_counts().head(10).index
    df_top = df[df['domain'].isin(top_domains)]
    
    pivot = df_top.pivot_table(values='impact_score', index='domain', 
                                columns='protein', aggfunc='mean')
    
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdYlGn_r', center=0.5,
                ax=ax7, cbar_kws={'label': 'Mean Impact Score'}, linewidths=0.5)
    ax7.set_title('Mean Impact Score by Protein & Domain (Top 10 Domains)', 
                  fontsize=12, fontweight='bold')
    ax7.set_xlabel('Protein', fontsize=11, fontweight='bold')
    ax7.set_ylabel('Domain', fontsize=11, fontweight='bold')
    
    # 8. Variant counts
    ax8 = fig.add_subplot(gs[2, 2])
    
    protein_counts = df['protein'].value_counts().sort_values(ascending=True)
    bars = ax8.barh(range(len(protein_counts)), protein_counts.values, 
                    color=colors_protein, edgecolor='black', alpha=0.7)
    ax8.set_yticks(range(len(protein_counts)))
    ax8.set_yticklabels(protein_counts.index)
    ax8.set_xlabel('Variant Count', fontsize=11, fontweight='bold')
    ax8.set_title('Variants per Protein', fontsize=12, fontweight='bold')
    ax8.grid(alpha=0.3, axis='x')
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, protein_counts.values)):
        ax8.text(val + 10, bar.get_y() + bar.get_height()/2, f'{val}',
                ha='left', va='center', fontsize=10, fontweight='bold')
    
    # Save figure
    output_file = '../plots/complete_dataset_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Visualization saved to: {output_file}")
    
    plt.close()


if __name__ == '__main__':
    csv_file = Path('../protein_csvs/domain_variants_complete.csv')

    if not csv_file.exists():
        print(f"Error: {csv_file} not found!")
        exit(1)

    visualize_complete_dataset(csv_file)
    print("\n✅ All done!")

