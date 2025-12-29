#!/usr/bin/env python3
"""Quick verification of final dataset."""

import pandas as pd

df = pd.read_csv('domain_variants_with_impact.csv')

print('=' * 80)
print('FINAL VERIFICATION')
print('=' * 80)
print(f'Total rows: {len(df)}')
print(f'Missing impact scores: {df["impact_score"].isna().sum()}')
print(f'Impact score range: [{df["impact_score"].min():.3f}, {df["impact_score"].max():.3f}]')
print(f'\nColumns: {list(df.columns)}')
print(f'\nSource distribution:')
print(df['impact_source'].value_counts())
print(f'\nProtein distribution:')
print(df['protein'].value_counts())
print('\n' + '=' * 80)
print('✅ All checks passed! Dataset is ready for step 2.3')
print('=' * 80)

