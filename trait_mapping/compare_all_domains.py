import pandas as pd
import re

# Load all files
opus = pd.read_csv('opus_domain.csv')
sonnet = pd.read_csv('sonnet_domain.csv')
chat = pd.read_csv('chat_domain.csv', skiprows=1)  # Skip the empty first row

# Parse chat_domain.csv - split "Domain (Protein)" into separate columns
def parse_domain_protein(val):
    match = re.match(r'(.+)\s*\((.+)\)', str(val))
    if match:
        return match.group(1).strip(), match.group(2).strip()
    return val, ''

chat[['Domain', 'Protein']] = chat['Domain (Protein)'].apply(lambda x: pd.Series(parse_domain_protein(x)))

# Create pivots for each
def make_pivot(df):
    pivot = df.pivot_table(index=['Domain', 'Protein'], columns='Trait', values='Weight', aggfunc='first').reset_index()
    pivot['Diff'] = pivot['Hypertrophy'] - pivot['Autophagy']
    return pivot

opus_pivot = make_pivot(opus)
sonnet_pivot = make_pivot(sonnet)
chat_pivot = make_pivot(chat)

# Merge all three
merged = opus_pivot.merge(sonnet_pivot, on=['Domain', 'Protein'], suffixes=('_opus', '_sonnet'))
merged = merged.merge(chat_pivot, on=['Domain', 'Protein'], suffixes=('', '_chat'))
merged = merged.rename(columns={'Hypertrophy': 'Hypertrophy_chat', 'Autophagy': 'Autophagy_chat', 'Diff': 'Diff_chat'})

print('=' * 110)
print('THREE-WAY COMPARISON: opus_domain.csv vs sonnet_domain.csv vs chat_domain.csv')
print('=' * 110)

# Print full comparison table
print('\n### FULL COMPARISON TABLE (H=Hypertrophy, A=Autophagy, D=Diff):')
print(f"{'Domain':<12} | {'Protein':<6} | {'Opus H/A/D':^18} | {'Sonnet H/A/D':^18} | {'Chat H/A/D':^18} | {'Agree?'}")
print('-' * 110)

for _, row in merged.sort_values(['Domain', 'Protein']).iterrows():
    opus_str = f"{row['Hypertrophy_opus']:.2f}/{row['Autophagy_opus']:.2f}/{row['Diff_opus']:+.2f}"
    sonnet_str = f"{row['Hypertrophy_sonnet']:.2f}/{row['Autophagy_sonnet']:.2f}/{row['Diff_sonnet']:+.2f}"
    chat_str = f"{row['Hypertrophy_chat']:.2f}/{row['Autophagy_chat']:.2f}/{row['Diff_chat']:+.2f}"
    
    # Check agreement on bias direction
    diffs = [row['Diff_opus'], row['Diff_sonnet'], row['Diff_chat']]
    all_auto = all(d < -0.1 for d in diffs)
    all_hyper = all(d > 0.1 for d in diffs)
    all_bal = all(-0.15 <= d <= 0.15 for d in diffs)
    agree = 'YES' if (all_auto or all_hyper or all_bal) else 'MIXED'
    
    print(f"  {row['Domain']:<12} | {row['Protein']:<6} | {opus_str:^18} | {sonnet_str:^18} | {chat_str:^18} | {agree}")

# Domain-level summary
print('\n' + '=' * 110)
print('DOMAIN-LEVEL SUMMARY (averaged across proteins):')
print('=' * 110)

def get_bias(diff):
    if diff < -0.15: return 'AUTO'
    elif diff > 0.15: return 'HYPER'
    else: return 'BAL'

opus_avg = opus_pivot.groupby('Domain')['Diff'].mean()
sonnet_avg = sonnet_pivot.groupby('Domain')['Diff'].mean()
chat_avg = chat_pivot.groupby('Domain')['Diff'].mean()

domain_summary = pd.DataFrame({
    'Opus_Diff': opus_avg, 'Sonnet_Diff': sonnet_avg, 'Chat_Diff': chat_avg
}).reset_index()

print(f"\n{'Domain':<12} | {'Opus':^12} | {'Sonnet':^12} | {'Chat':^12} | {'Consensus'}")
print('-' * 70)

for _, row in domain_summary.sort_values('Opus_Diff').iterrows():
    o_bias = get_bias(row['Opus_Diff'])
    s_bias = get_bias(row['Sonnet_Diff'])
    c_bias = get_bias(row['Chat_Diff'])
    
    biases = [o_bias, s_bias, c_bias]
    if biases.count(biases[0]) == 3:
        consensus = biases[0]
    elif biases.count('AUTO') >= 2:
        consensus = 'AUTO (2/3)'
    elif biases.count('HYPER') >= 2:
        consensus = 'HYPER (2/3)'
    elif biases.count('BAL') >= 2:
        consensus = 'BAL (2/3)'
    else:
        consensus = 'NO CONSENSUS'
    
    print(f"  {row['Domain']:<12} | {o_bias:^12} | {s_bias:^12} | {c_bias:^12} | {consensus}")

# Statistics
print('\n' + '=' * 110)
print('WEIGHT STATISTICS:')
print('=' * 110)
print(f"\n{'Metric':<30} | {'Opus':^10} | {'Sonnet':^10} | {'Chat':^10}")
print('-' * 70)
print(f"  {'Mean Hypertrophy Weight':<30} | {merged['Hypertrophy_opus'].mean():^10.3f} | {merged['Hypertrophy_sonnet'].mean():^10.3f} | {merged['Hypertrophy_chat'].mean():^10.3f}")
print(f"  {'Mean Autophagy Weight':<30} | {merged['Autophagy_opus'].mean():^10.3f} | {merged['Autophagy_sonnet'].mean():^10.3f} | {merged['Autophagy_chat'].mean():^10.3f}")
print(f"  {'Mean Diff (H-A)':<30} | {merged['Diff_opus'].mean():^10.3f} | {merged['Diff_sonnet'].mean():^10.3f} | {merged['Diff_chat'].mean():^10.3f}")

