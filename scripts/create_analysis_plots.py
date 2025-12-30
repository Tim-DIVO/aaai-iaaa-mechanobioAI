#!/usr/bin/env python3
"""
Create analysis plots for:
1. Variant impact scores (observed vs imputed)
2. Variant allele frequencies (observed vs imputed)
3. Synthetic cohort distributions
4. Model training and evaluation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150

# Load data
variants = pd.read_csv('../protein_csvs/domain_variants_complete.csv')
cohort = pd.read_csv('../synthetic_cohort/synthetic_cohort_summary.csv')
phenotypes = pd.read_csv('../trait_mapping/phenotype_definitions.csv')

print("Creating analysis plots...")

# =============================================================================
# PLOT 1: Impact Score Distribution (Observed vs Imputed)
# =============================================================================
fig1, axes = plt.subplots(1, 2, figsize=(12, 5))

# Split by source
observed = variants[variants['impact_source'] == 'both']
imputed = variants[variants['impact_source'] == 'sampled']
sift_only = variants[variants['impact_source'] == 'sift_only'] if 'sift_only' in variants['impact_source'].values else pd.DataFrame()

ax1 = axes[0]
ax1.hist(observed['impact_score'], bins=30, alpha=0.7, label=f'Observed (n={len(observed)})', color='#2E7D32', density=True)
ax1.hist(imputed['impact_score'], bins=30, alpha=0.7, label=f'Imputed (n={len(imputed)})', color='#E65100', density=True)
ax1.set_xlabel('Impact Score')
ax1.set_ylabel('Density')
ax1.set_title('Impact Score Distribution\n(Observed vs Imputed)')
ax1.legend()

# By consequence type for imputed
ax2 = axes[1]
consequence_order = ['missense', 'stop gained', 'frameshift', 'inframe deletion']
colors = {'missense': '#1976D2', 'stop gained': '#D32F2F', 'frameshift': '#7B1FA2', 'inframe deletion': '#F57C00'}
for cons in consequence_order:
    subset = imputed[imputed['consequence_type'] == cons]
    if len(subset) > 0:
        ax2.hist(subset['impact_score'], bins=20, alpha=0.6, label=f'{cons} (n={len(subset)})', 
                 color=colors.get(cons, 'gray'), density=True)
ax2.set_xlabel('Impact Score')
ax2.set_ylabel('Density')
ax2.set_title('Imputed Impact Scores\nby Consequence Type')
ax2.legend(fontsize=8)

fig1.tight_layout()
fig1.savefig('../plots/variant_impact_scores.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved variant_impact_scores.png")
plt.close(fig1)

# =============================================================================
# PLOT 2: Allele Frequency Distribution (Observed vs Imputed)
# =============================================================================
fig2, axes = plt.subplots(1, 2, figsize=(12, 5))

obs_freq = variants[variants['frequency_source'] == 'observed']
imp_freq = variants[variants['frequency_source'] == 'sampled']

# Filter out zero frequencies for log transform
obs_freq_nonzero = obs_freq[obs_freq['allele_frequency'] > 0]
imp_freq_nonzero = imp_freq[imp_freq['allele_frequency'] > 0]

ax1 = axes[0]
ax1.hist(np.log10(obs_freq_nonzero['allele_frequency']), bins=30, alpha=0.7,
         label=f'Observed (n={len(obs_freq_nonzero)})', color='#2E7D32', density=True)
ax1.hist(np.log10(imp_freq_nonzero['allele_frequency']), bins=30, alpha=0.7,
         label=f'Imputed (n={len(imp_freq_nonzero)})', color='#E65100', density=True)
ax1.set_xlabel('log₁₀(Allele Frequency)')
ax1.set_ylabel('Density')
ax1.set_title('Allele Frequency Distribution\n(Observed vs Imputed)')
ax1.legend()

# By protein
ax2 = axes[1]
protein_colors = {'TTN': '#1976D2', 'MuRF1': '#D32F2F', 'MuRF2': '#7B1FA2', 'p62': '#388E3C', 'NBR1': '#F57C00'}
for protein in ['TTN', 'MuRF1', 'MuRF2', 'p62', 'NBR1']:
    subset = variants[(variants['protein'] == protein) & (variants['allele_frequency'] > 0)]
    ax2.hist(np.log10(subset['allele_frequency']), bins=25, alpha=0.5,
             label=f'{protein} (n={len(subset)})', color=protein_colors[protein], density=True)
ax2.set_xlabel('log₁₀(Allele Frequency)')
ax2.set_ylabel('Density')
ax2.set_title('Allele Frequency by Protein')
ax2.legend(fontsize=8)

fig2.tight_layout()
fig2.savefig('../plots/variant_allele_frequencies.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved variant_allele_frequencies.png")
plt.close(fig2)

# =============================================================================
# PLOT 3: Synthetic Cohort Phenotype Distributions
# =============================================================================
fig3, axes = plt.subplots(2, 5, figsize=(16, 6))
pheno_cols = phenotypes['phenotype'].tolist()

for i, (ax, pheno) in enumerate(zip(axes.flat, pheno_cols)):
    ax.hist(cohort[pheno], bins=50, color=plt.cm.tab10(i), alpha=0.7, edgecolor='white')
    ax.axvline(cohort[pheno].mean(), color='black', linestyle='--', linewidth=1.5, label=f'μ={cohort[pheno].mean():.2f}')
    ax.set_xlabel('Value (z-score)')
    ax.set_title(pheno.replace('_', ' ').title(), fontsize=10)
    ax.legend(fontsize=7)

fig3.suptitle('Synthetic Cohort: 10 Phenotype Distributions (n=50,000)', fontsize=14, fontweight='bold')
fig3.tight_layout()
fig3.savefig('../plots/cohort_phenotype_distributions.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved cohort_phenotype_distributions.png")
plt.close(fig3)

# =============================================================================
# PLOT 4: Cohort Burden Score Scatter with Phenotype Coloring
# =============================================================================
fig4, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# Sample for visibility
sample_idx = np.random.choice(len(cohort), size=5000, replace=False)
sample = cohort.iloc[sample_idx]

ax1 = axes[0]
scatter = ax1.scatter(sample['hypertrophy_burden'], sample['autophagy_burden'], 
                      c=sample['natural_muscle_mass'], cmap='RdYlBu', alpha=0.5, s=10)
ax1.set_xlabel('Hypertrophy Burden')
ax1.set_ylabel('Autophagy Burden')
ax1.set_title('Burden Scores colored by\nNatural Muscle Mass')
plt.colorbar(scatter, ax=ax1, label='Muscle Mass')

ax2 = axes[1]
scatter2 = ax2.scatter(sample['hypertrophy_burden'], sample['autophagy_burden'],
                       c=sample['recovery_speed'], cmap='RdYlBu', alpha=0.5, s=10)
ax2.set_xlabel('Hypertrophy Burden')
ax2.set_ylabel('Autophagy Burden')
ax2.set_title('Burden Scores colored by\nRecovery Speed')
plt.colorbar(scatter2, ax=ax2, label='Recovery Speed')

ax3 = axes[2]
scatter3 = ax3.scatter(sample['hypertrophy_burden'], sample['autophagy_burden'],
                       c=sample['training_responsiveness'], cmap='RdYlBu', alpha=0.5, s=10)
ax3.set_xlabel('Hypertrophy Burden')
ax3.set_ylabel('Autophagy Burden')
ax3.set_title('Burden Scores colored by\nTraining Responsiveness')
plt.colorbar(scatter3, ax=ax3, label='Training Resp.')

fig4.suptitle('Burden Score Space with Phenotype Coloring', fontsize=12, fontweight='bold')
fig4.tight_layout()
fig4.savefig('../plots/cohort_burden_phenotype_scatter.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved cohort_burden_phenotype_scatter.png")
plt.close(fig4)

print("\nVariant and cohort plots complete!")

