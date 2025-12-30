#!/usr/bin/env python3
"""
Create 5 publication-worthy figures for conference poster.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches

# Style settings
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['figure.dpi'] = 150

# Load data
cohort = pd.read_csv('../synthetic_cohort/synthetic_cohort_summary.csv')
variants = pd.read_csv('../protein_csvs/domain_variants_complete.csv')
weights = pd.read_csv('../trait_mapping/combined_domain.csv')
phenotypes = pd.read_csv('../trait_mapping/phenotype_definitions.csv')

print("Creating 5 poster figures...")

# =============================================================================
# FIGURE 1: Pipeline Overview / Knowledge Graph Structure
# =============================================================================
fig1, ax1 = plt.subplots(figsize=(12, 6))
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 5)
ax1.axis('off')
ax1.set_title('MechanobioAI: Genotype-to-Phenotype Pipeline', fontsize=16, fontweight='bold', pad=20)

# Boxes
boxes = [
    (0.5, 2, 1.8, 1.5, '#E8F4FD', 'Variants\n(n=2,139)', 'SNPs in titin\nsignalosome'),
    (3, 2, 1.8, 1.5, '#FFF3E0', 'Domains\n(n=21)', 'Structured +\nIDR regions'),
    (5.5, 2, 1.8, 1.5, '#E8F5E9', 'Pathway\nBurdens', 'Hypertrophy\nAutophagy'),
    (8, 2, 1.8, 1.5, '#FCE4EC', 'Phenotypes\n(n=10)', 'Observable\ntraits'),
]

for x, y, w, h, color, title, subtitle in boxes:
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", 
                          facecolor=color, edgecolor='#333', linewidth=2)
    ax1.add_patch(rect)
    ax1.text(x + w/2, y + h*0.7, title, ha='center', va='center', fontsize=11, fontweight='bold')
    ax1.text(x + w/2, y + h*0.3, subtitle, ha='center', va='center', fontsize=9, color='#555')

# Arrows
arrow_props = dict(arrowstyle='->', color='#333', lw=2)
ax1.annotate('', xy=(3, 2.75), xytext=(2.3, 2.75), arrowprops=arrow_props)
ax1.annotate('', xy=(5.5, 2.75), xytext=(4.8, 2.75), arrowprops=arrow_props)
ax1.annotate('', xy=(8, 2.75), xytext=(7.3, 2.75), arrowprops=arrow_props)

# Labels on arrows
ax1.text(2.65, 3.1, 'impact\nscore', ha='center', fontsize=8, color='#555')
ax1.text(5.15, 3.1, 'domain\nweights', ha='center', fontsize=8, color='#555')
ax1.text(7.65, 3.1, 'pathway\ncoefficients', ha='center', fontsize=8, color='#555')

# Proteins at bottom
proteins = ['TTN', 'MuRF1', 'MuRF2', 'p62', 'NBR1']
for i, p in enumerate(proteins):
    ax1.text(0.5 + i*0.4, 1.3, p, fontsize=8, ha='center', color='#1565C0')

ax1.text(1.4, 0.8, '5 proteins in titin kinase signalosome', fontsize=9, ha='center', style='italic')

fig1.tight_layout()
fig1.savefig('../plots/fig1_pipeline_overview.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved fig1_pipeline_overview.png")
plt.close(fig1)

# =============================================================================
# FIGURE 2: Domain-Trait Weight Heatmap
# =============================================================================
fig2, ax2 = plt.subplots(figsize=(10, 8))

# Pivot weights for heatmap
weights_hyp = weights[weights['Trait'] == 'Hypertrophy'][['Domain', 'Protein', 'Weight']]
weights_auto = weights[weights['Trait'] == 'Autophagy'][['Domain', 'Protein', 'Weight']]
weights_hyp['domain_protein'] = weights_hyp['Domain'] + ' (' + weights_hyp['Protein'] + ')'
weights_auto['domain_protein'] = weights_auto['Domain'] + ' (' + weights_auto['Protein'] + ')'

heatmap_data = pd.DataFrame({
    'Domain': weights_hyp['domain_protein'].values,
    'Hypertrophy': weights_hyp['Weight'].values,
    'Autophagy': weights_auto['Weight'].values
}).set_index('Domain')

sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap='RdYlBu_r', 
            center=0.5, vmin=0, vmax=1, ax=ax2,
            cbar_kws={'label': 'Pathway Involvement Weight'})
ax2.set_title('Domain-Pathway Weight Matrix\n(LLM-assisted literature review)', fontsize=14, fontweight='bold')
ax2.set_xlabel('')
ax2.set_ylabel('')
plt.xticks(rotation=0)
plt.yticks(rotation=0)

fig2.tight_layout()
fig2.savefig('../plots/fig2_domain_weights_heatmap.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved fig2_domain_weights_heatmap.png")
plt.close(fig2)

# =============================================================================
# FIGURE 3: Phenotype-Pathway Coefficient Plot (showing differentiation)
# =============================================================================
fig3, ax3 = plt.subplots(figsize=(10, 6))

pheno_names = phenotypes['phenotype'].tolist()
hyp_coefs = phenotypes['hypertrophy_coef'].values
auto_coefs = phenotypes['autophagy_coef'].values

x = np.arange(len(pheno_names))
width = 0.35

bars1 = ax3.bar(x - width/2, hyp_coefs, width, label='Hypertrophy Pathway', color='#E53935', alpha=0.8)
bars2 = ax3.bar(x + width/2, auto_coefs, width, label='Autophagy Pathway', color='#1E88E5', alpha=0.8)

ax3.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax3.set_ylabel('Pathway Coefficient', fontsize=12)
ax3.set_title('Phenotype-Pathway Coefficients\n(Disruption effects on observable traits)', fontsize=14, fontweight='bold')
ax3.set_xticks(x)
ax3.set_xticklabels([p.replace('_', '\n') for p in pheno_names], rotation=45, ha='right', fontsize=9)
ax3.legend(loc='upper right')

# Add interpretation annotations
ax3.annotate('Muscle phenotypes:\nHypertrophy-dominant', xy=(1.5, -1.3), fontsize=9,
             ha='center', style='italic', color='#E53935')
ax3.annotate('Endurance phenotypes:\nAutophagy-dominant', xy=(6, -1.3), fontsize=9,
             ha='center', style='italic', color='#1E88E5')

fig3.tight_layout()
fig3.savefig('../plots/fig3_phenotype_coefficients.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved fig3_phenotype_coefficients.png")
plt.close(fig3)

# =============================================================================
# FIGURE 4: Partial Correlation Plot (unique pathway contributions)
# =============================================================================
def partial_corr(x, y, z):
    '''Partial correlation of x and y controlling for z'''
    slope_xz = np.cov(x, z)[0,1] / np.var(z)
    slope_yz = np.cov(y, z)[0,1] / np.var(z)
    x_resid = x - slope_xz * z
    y_resid = y - slope_yz * z
    return np.corrcoef(x_resid, y_resid)[0,1]

hyp = cohort['hypertrophy_burden'].values
auto = cohort['autophagy_burden'].values

partial_hyp = []
partial_auto = []
for pheno in pheno_names:
    p = cohort[pheno].values
    partial_hyp.append(partial_corr(hyp, p, auto))
    partial_auto.append(partial_corr(auto, p, hyp))

fig4, ax4 = plt.subplots(figsize=(10, 8))

# Scatter plot of partial correlations
colors = ['#E53935' if abs(h) > abs(a) + 0.1 else '#1E88E5' if abs(a) > abs(h) + 0.1 else '#9C27B0'
          for h, a in zip(partial_hyp, partial_auto)]

for i, (h, a, name, c) in enumerate(zip(partial_hyp, partial_auto, pheno_names, colors)):
    ax4.scatter(h, a, s=200, c=c, alpha=0.8, edgecolors='black', linewidths=1)
    ax4.annotate(name.replace('_', ' '), (h, a), xytext=(5, 5), textcoords='offset points', fontsize=9)

ax4.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
ax4.axvline(x=0, color='gray', linestyle='--', linewidth=0.5)
ax4.plot([-1, 1], [-1, 1], 'k:', alpha=0.3, label='Equal contribution')

ax4.set_xlabel('Partial r (Hypertrophy | Autophagy)', fontsize=12)
ax4.set_ylabel('Partial r (Autophagy | Hypertrophy)', fontsize=12)
ax4.set_title('Unique Pathway Contributions to Phenotypes\n(Partial correlations controlling for the other pathway)',
              fontsize=14, fontweight='bold')

# Legend
hyp_patch = mpatches.Patch(color='#E53935', label='Hypertrophy-dominant')
auto_patch = mpatches.Patch(color='#1E88E5', label='Autophagy-dominant')
both_patch = mpatches.Patch(color='#9C27B0', label='Both pathways')
ax4.legend(handles=[hyp_patch, auto_patch, both_patch], loc='lower left')

ax4.set_xlim(-1.1, 1.1)
ax4.set_ylim(-1.1, 1.1)

fig4.tight_layout()
fig4.savefig('../plots/fig4_partial_correlations.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved fig4_partial_correlations.png")
plt.close(fig4)

# =============================================================================
# FIGURE 5: Burden Distribution & Phenotype Scatter
# =============================================================================
fig5, axes = plt.subplots(1, 3, figsize=(14, 4.5))

# Panel A: Burden distributions
ax5a = axes[0]
ax5a.hist(cohort['hypertrophy_burden'], bins=50, alpha=0.7, color='#E53935', label='Hypertrophy', density=True)
ax5a.hist(cohort['autophagy_burden'], bins=50, alpha=0.7, color='#1E88E5', label='Autophagy', density=True)
ax5a.set_xlabel('Pathway Burden Score')
ax5a.set_ylabel('Density')
ax5a.set_title('A. Burden Score Distributions\n(n=50,000 synthetic individuals)')
ax5a.legend()

# Panel B: Muscle phenotype vs hypertrophy (with autophagy color)
ax5b = axes[1]
scatter = ax5b.scatter(cohort['hypertrophy_burden'], cohort['natural_muscle_mass'],
                       c=cohort['autophagy_burden'], cmap='coolwarm', alpha=0.3, s=1)
ax5b.set_xlabel('Hypertrophy Burden')
ax5b.set_ylabel('Natural Muscle Mass (z-score)')
ax5b.set_title('B. Muscle Mass vs Hypertrophy Burden\n(colored by Autophagy burden)')
plt.colorbar(scatter, ax=ax5b, label='Autophagy Burden')

# Panel C: Recovery phenotype vs autophagy (with hypertrophy color)
ax5c = axes[2]
scatter2 = ax5c.scatter(cohort['autophagy_burden'], cohort['recovery_speed'],
                        c=cohort['hypertrophy_burden'], cmap='coolwarm', alpha=0.3, s=1)
ax5c.set_xlabel('Autophagy Burden')
ax5c.set_ylabel('Recovery Speed (z-score)')
ax5c.set_title('C. Recovery vs Autophagy Burden\n(colored by Hypertrophy burden)')
plt.colorbar(scatter2, ax=ax5c, label='Hypertrophy Burden')

fig5.tight_layout()
fig5.savefig('../plots/fig5_burden_phenotype_scatter.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved fig5_burden_phenotype_scatter.png")
plt.close(fig5)

print("\nAll 5 figures saved to ../plots/")

