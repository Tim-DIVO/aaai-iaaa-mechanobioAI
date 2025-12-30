#!/usr/bin/env python3
"""
Generate synthetic cohort with genotypes, disruption scores, and 10 observable phenotypes.

Burden scores:
    hypertrophy_burden = Σ (impact_score × hypertrophy_weight × genotype_dosage)
    autophagy_burden = Σ (impact_score × autophagy_weight × genotype_dosage)

Phenotypes are linear combinations of burden scores with trait-specific coefficients:
    phenotype = baseline + hyp_coef * hyp_burden + auto_coef * auto_burden + noise

Phenotype mappings based on exercise physiology literature.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Config
N_INDIVIDUALS = 50000
SEED = 42
OUTPUT_DIR = '../synthetic_cohort'

np.random.seed(SEED)

# Load data
variants = pd.read_csv('../protein_csvs/domain_variants_complete.csv')
weights = pd.read_csv('../trait_mapping/combined_domain.csv')

# Create weight lookup: (domain, protein, trait) -> weight
weight_lookup = {}
for _, row in weights.iterrows():
    key = (row['Domain'], row['Protein'], row['Trait'])
    weight_lookup[key] = row['Weight']

# Get hypertrophy and autophagy weights for each variant
def get_weight(row, trait):
    key = (row['domain'], row['protein'], trait)
    return weight_lookup.get(key, 0.0)

variants['hypertrophy_weight'] = variants.apply(lambda r: get_weight(r, 'Hypertrophy'), axis=1)
variants['autophagy_weight'] = variants.apply(lambda r: get_weight(r, 'Autophagy'), axis=1)

# Create variant info for output
variant_ids = [f"{r['protein']}_{r['domain']}_{r['position']}_{r['wild_type']}{r['mutated_type']}" 
               for _, r in variants.iterrows()]
variants['variant_id'] = variant_ids

print(f"Loaded {len(variants)} variants")
print(f"Generating {N_INDIVIDUALS} synthetic individuals...")

# Generate genotype matrix: N_INDIVIDUALS x N_VARIANTS
# Each genotype is 0, 1, or 2 based on allele frequency
allele_freqs = variants['allele_frequency'].values
genotype_matrix = np.random.binomial(2, allele_freqs, size=(N_INDIVIDUALS, len(variants)))

print(f"Genotype matrix shape: {genotype_matrix.shape}")
print(f"Non-zero genotypes: {np.sum(genotype_matrix > 0)} ({100*np.mean(genotype_matrix > 0):.2f}%)")

# Compute burden scores for each individual
impact_scores = variants['impact_score'].values
hyp_weights = variants['hypertrophy_weight'].values
auto_weights = variants['autophagy_weight'].values

# burden = genotype × impact × weight (summed over variants)
hypertrophy_burden = genotype_matrix @ (impact_scores * hyp_weights)
autophagy_burden = genotype_matrix @ (impact_scores * auto_weights)

print(f"\nBurden statistics:")
print(f"  Hypertrophy burden: mean={hypertrophy_burden.mean():.4f}, std={hypertrophy_burden.std():.4f}")
print(f"  Autophagy burden: mean={autophagy_burden.mean():.4f}, std={autophagy_burden.std():.4f}")

# Load phenotype definitions
phenotypes = pd.read_csv('../trait_mapping/phenotype_definitions.csv')
print(f"\nGenerating {len(phenotypes)} phenotypes...")

# Create output dataframe with burden scores
summary_df = pd.DataFrame({
    'individual_id': range(N_INDIVIDUALS),
    'hypertrophy_burden': hypertrophy_burden,
    'autophagy_burden': autophagy_burden,
    'n_variants': np.sum(genotype_matrix > 0, axis=1),
    'total_dosage': np.sum(genotype_matrix, axis=1)
})

# Generate each phenotype
for _, pheno in phenotypes.iterrows():
    name = pheno['phenotype']
    hyp_coef = pheno['hypertrophy_coef']
    auto_coef = pheno['autophagy_coef']
    noise_std = pheno['noise_std']
    baseline = pheno['baseline']

    # phenotype = baseline + hyp_coef * hyp_burden + auto_coef * auto_burden + noise
    noise = np.random.normal(0, noise_std, N_INDIVIDUALS)
    values = baseline + hyp_coef * hypertrophy_burden + auto_coef * autophagy_burden + noise
    summary_df[name] = values

    print(f"  {name}: mean={values.mean():.3f}, std={values.std():.3f}")

# Save summary
summary_df.to_csv(f'{OUTPUT_DIR}/synthetic_cohort_summary.csv', index=False)
print(f"\nSaved summary to {OUTPUT_DIR}/synthetic_cohort_summary.csv")

# Save full genotype matrix (sparse representation for efficiency)
print("Saving genotype records (this may take a moment)...")
genotype_records = []
for i in range(N_INDIVIDUALS):
    for j in np.where(genotype_matrix[i] > 0)[0]:
        genotype_records.append({
            'individual_id': i,
            'variant_id': variants.iloc[j]['variant_id'],
            'genotype': genotype_matrix[i, j]
        })

genotype_df = pd.DataFrame(genotype_records)
genotype_df.to_csv(f'{OUTPUT_DIR}/synthetic_genotypes.csv', index=False)
print(f"Saved {len(genotype_df)} genotype records to {OUTPUT_DIR}/synthetic_genotypes.csv")

# Print correlations
print(f"\n=== Phenotype Correlations with Burden Scores ===")
pheno_cols = phenotypes['phenotype'].tolist()
for col in pheno_cols:
    r_hyp = np.corrcoef(hypertrophy_burden, summary_df[col])[0,1]
    r_auto = np.corrcoef(autophagy_burden, summary_df[col])[0,1]
    print(f"  {col}: r_hyp={r_hyp:.3f}, r_auto={r_auto:.3f}")

