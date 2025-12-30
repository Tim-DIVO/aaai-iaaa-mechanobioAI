import pandas as pd

# Load both files
opus = pd.read_csv('opus_domain.csv')
claude = pd.read_csv('claude_domain.csv')

# Create pivots
def make_pivot(df):
    pivot = df.pivot_table(index=['Domain', 'Protein'], columns='Trait', values='Weight', aggfunc='first').reset_index()
    pivot['Diff'] = pivot['Hypertrophy'] - pivot['Autophagy']
    return pivot

opus_pivot = make_pivot(opus)
claude_pivot = make_pivot(claude)

# Merge for comparison
merged = opus_pivot.merge(claude_pivot, on=['Domain', 'Protein'], suffixes=('_opus', '_claude'))

print('=' * 100)
print('COMPARISON: opus_domain.csv vs claude_domain.csv')
print('=' * 100)

# Calculate agreement/disagreement
merged['H_diff'] = merged['Hypertrophy_opus'] - merged['Hypertrophy_claude']
merged['A_diff'] = merged['Autophagy_opus'] - merged['Autophagy_claude']
merged['Bias_diff'] = merged['Diff_opus'] - merged['Diff_claude']  # Difference in H-A bias

print('\n### SIMILAR ASSESSMENTS (Both Hypertrophy and Autophagy weights within 0.15):')
similar = merged[(abs(merged['H_diff']) <= 0.15) & (abs(merged['A_diff']) <= 0.15)]
similar = similar.sort_values('Domain')
for _, row in similar.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | "
          f"Opus: H={row['Hypertrophy_opus']:.2f}/A={row['Autophagy_opus']:.2f} | "
          f"Claude: H={row['Hypertrophy_claude']:.2f}/A={row['Autophagy_claude']:.2f}")

print('\n### HYPERTROPHY WEIGHT DIFFERENCES (Opus higher by >0.15):')
h_opus_higher = merged[merged['H_diff'] > 0.15].sort_values('H_diff', ascending=False)
for _, row in h_opus_higher.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | "
          f"Opus H={row['Hypertrophy_opus']:.2f} vs Claude H={row['Hypertrophy_claude']:.2f} | Delta={row['H_diff']:+.2f}")

print('\n### HYPERTROPHY WEIGHT DIFFERENCES (Claude higher by >0.15):')
h_claude_higher = merged[merged['H_diff'] < -0.15].sort_values('H_diff')
for _, row in h_claude_higher.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | "
          f"Opus H={row['Hypertrophy_opus']:.2f} vs Claude H={row['Hypertrophy_claude']:.2f} | Delta={row['H_diff']:+.2f}")

print('\n### AUTOPHAGY WEIGHT DIFFERENCES (Opus higher by >0.15):')
a_opus_higher = merged[merged['A_diff'] > 0.15].sort_values('A_diff', ascending=False)
for _, row in a_opus_higher.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | "
          f"Opus A={row['Autophagy_opus']:.2f} vs Claude A={row['Autophagy_claude']:.2f} | Delta={row['A_diff']:+.2f}")

print('\n### AUTOPHAGY WEIGHT DIFFERENCES (Claude higher by >0.15):')
a_claude_higher = merged[merged['A_diff'] < -0.15].sort_values('A_diff')
for _, row in a_claude_higher.iterrows():
    print(f"  {row['Domain']:12} | {row['Protein']:6} | "
          f"Opus A={row['Autophagy_opus']:.2f} vs Claude A={row['Autophagy_claude']:.2f} | Delta={row['A_diff']:+.2f}")

print('\n' + '=' * 100)
print('DOMAIN-LEVEL SUMMARY (averaged across proteins):')
print('=' * 100)

opus_avg = opus_pivot.groupby('Domain').agg({'Hypertrophy': 'mean', 'Autophagy': 'mean', 'Diff': 'mean'}).reset_index()
claude_avg = claude_pivot.groupby('Domain').agg({'Hypertrophy': 'mean', 'Autophagy': 'mean', 'Diff': 'mean'}).reset_index()
domain_compare = opus_avg.merge(claude_avg, on='Domain', suffixes=('_opus', '_claude'))
domain_compare['Bias_agreement'] = (domain_compare['Diff_opus'] * domain_compare['Diff_claude']) > 0  # Same sign = agree

print('\nDomain         | Opus Bias  | Claude Bias | Agreement')
print('-' * 60)
for _, row in domain_compare.sort_values('Diff_opus').iterrows():
    opus_bias = 'AUTO' if row['Diff_opus'] < -0.15 else ('HYPER' if row['Diff_opus'] > 0.15 else 'BAL')
    claude_bias = 'AUTO' if row['Diff_claude'] < -0.15 else ('HYPER' if row['Diff_claude'] > 0.15 else 'BAL')
    agree = 'YES' if opus_bias == claude_bias else 'DIFFERS'
    print(f"  {row['Domain']:12} | {opus_bias:10} | {claude_bias:11} | {agree}")

