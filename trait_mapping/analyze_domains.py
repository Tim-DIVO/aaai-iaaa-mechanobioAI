import pandas as pd
import sys

filename = sys.argv[1] if len(sys.argv) > 1 else 'opus_domain.csv'
df = pd.read_csv(filename)

# Pivot to compare Hypertrophy vs Autophagy weights for each Domain-Protein pair
pivot = df.pivot_table(index=['Domain', 'Protein'], columns='Trait', values='Weight', aggfunc='first').reset_index()
pivot['Diff'] = pivot['Hypertrophy'] - pivot['Autophagy']
pivot['Abs_Diff'] = abs(pivot['Diff'])

# Sort by absolute difference
pivot_sorted = pivot.sort_values('Abs_Diff', ascending=False)

print('=' * 80)
print(f'DIFFERENTIAL ANALYSIS: {filename}')
print('Diff > 0 means MORE important for Hypertrophy')
print('Diff < 0 means MORE important for Autophagy')
print('=' * 80)

print('\n### MOST DIFFERENTIAL (Autophagy-biased, Diff < -0.4):')
autophagy_biased = pivot_sorted[pivot_sorted['Diff'] < -0.4]
for _, row in autophagy_biased.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | H={row['Hypertrophy']:.2f} | A={row['Autophagy']:.2f} | Diff={row['Diff']:+.2f}")

print('\n### MODERATELY DIFFERENTIAL (Autophagy-biased, -0.4 < Diff < -0.15):')
mod_autophagy = pivot_sorted[(pivot_sorted['Diff'] >= -0.4) & (pivot_sorted['Diff'] < -0.15)]
for _, row in mod_autophagy.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | H={row['Hypertrophy']:.2f} | A={row['Autophagy']:.2f} | Diff={row['Diff']:+.2f}")

print('\n### FAIRLY EQUAL (-0.15 <= Diff <= +0.15):')
equal = pivot_sorted[(pivot_sorted['Diff'] >= -0.15) & (pivot_sorted['Diff'] <= 0.15)]
for _, row in equal.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | H={row['Hypertrophy']:.2f} | A={row['Autophagy']:.2f} | Diff={row['Diff']:+.2f}")

print('\n### MODERATELY DIFFERENTIAL (Hypertrophy-biased, +0.15 < Diff):')
mod_hyper = pivot_sorted[pivot_sorted['Diff'] > 0.15]
for _, row in mod_hyper.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | H={row['Hypertrophy']:.2f} | A={row['Autophagy']:.2f} | Diff={row['Diff']:+.2f}")

print('\n' + '=' * 80)
print('SUMMARY BY DOMAIN (averaged across proteins):')
print('=' * 80)
domain_avg = pivot.groupby('Domain').agg({'Hypertrophy': 'mean', 'Autophagy': 'mean'}).reset_index()
domain_avg['Diff'] = domain_avg['Hypertrophy'] - domain_avg['Autophagy']
domain_avg = domain_avg.sort_values('Diff')
for _, row in domain_avg.iterrows():
    bias = 'AUTOPHAGY' if row['Diff'] < -0.2 else ('HYPERTROPHY' if row['Diff'] > 0.2 else 'BALANCED')
    print(f"  {row['Domain']:12} | H={row['Hypertrophy']:.2f} | A={row['Autophagy']:.2f} | Diff={row['Diff']:+.2f} | {bias}")

