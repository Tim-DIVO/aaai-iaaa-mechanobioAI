#!/usr/bin/env python3
"""
Filter variants CSV to only retain variants that are in a domain.
Removes all variants where domain == "None".
"""

import pandas as pd
import sys
from pathlib import Path


def filter_domain_variants(input_csv: Path, output_csv: Path):
    """Filter CSV to only keep variants that are in a domain."""
    
    print(f"Loading variants from {input_csv}...")
    df = pd.read_csv(input_csv)
    
    print(f"  Total rows: {len(df)}")
    print(f"  Unique proteins: {df['protein'].nunique()}")
    
    # Count variants before filtering
    total_variants = len(df)
    none_domain_count = df['domain'].isna().sum()

    print(f"\nFiltering out variants with no domain (NaN)...")
    print(f"  Variants with no domain: {none_domain_count}")

    # Filter to only keep variants in domains (remove NaN values)
    df_filtered = df[df['domain'].notna()].copy()
    
    print(f"  Variants in domains: {len(df_filtered)}")
    
    # Write filtered CSV
    print(f"\nWriting filtered variants to {output_csv}...")
    df_filtered.to_csv(output_csv, index=False)
    
    # Summary statistics
    print(f"\n{'=' * 80}")
    print("SUMMARY")
    print(f"{'=' * 80}")
    print(f"  Input rows: {total_variants}")
    print(f"  Output rows: {len(df_filtered)}")
    print(f"  Removed rows: {total_variants - len(df_filtered)}")
    print(f"  Retention rate: {100 * len(df_filtered) / total_variants:.1f}%")
    
    print(f"\nBreakdown by protein:")
    for protein in sorted(df['protein'].unique()):
        total = len(df[df['protein'] == protein])
        in_domain = len(df_filtered[df_filtered['protein'] == protein])
        print(f"  {protein:8s}: {in_domain:5d} / {total:5d} ({100 * in_domain / total:.1f}%)")
    
    print(f"\n✅ Filtered CSV saved to: {output_csv}")
    print(f"{'=' * 80}")


if __name__ == '__main__':
    input_file = Path('../protein_csvs/all_variants.csv')
    output_file = Path('../protein_csvs/domain_variants.csv')

    if not input_file.exists():
        print(f"Error: {input_file} not found!")
        sys.exit(1)

    filter_domain_variants(input_file, output_file)

