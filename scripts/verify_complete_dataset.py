#!/usr/bin/env python3
"""
Final verification of the complete dataset before proceeding to Step 2.3.
"""

import pandas as pd
import numpy as np
from pathlib import Path


def verify_complete_dataset(csv_file: Path):
    """Verify the complete dataset is ready for Step 2.3."""
    
    print("=" * 80)
    print("FINAL DATASET VERIFICATION")
    print("=" * 80)
    print(f"\nLoading: {csv_file}")
    
    df = pd.read_csv(csv_file)
    
    print(f"Total rows: {len(df)}")
    print(f"Columns: {list(df.columns)}")
    
    # Check for required columns
    required_cols = ['protein', 'domain', 'position', 'wild_type', 'mutated_type',
                     'consequence_type', 'impact_score', 'allele_frequency',
                     'impact_source', 'frequency_source']
    
    print("\n" + "=" * 80)
    print("REQUIRED COLUMNS CHECK")
    print("=" * 80)
    
    all_present = True
    for col in required_cols:
        present = col in df.columns
        status = "✅" if present else "❌"
        print(f"  {status} {col}")
        if not present:
            all_present = False
    
    if not all_present:
        print("\n❌ FAILED: Missing required columns!")
        return False
    
    # Check for missing values
    print("\n" + "=" * 80)
    print("MISSING VALUES CHECK")
    print("=" * 80)
    
    critical_cols = ['protein', 'domain', 'impact_score', 'allele_frequency']
    no_missing = True
    
    for col in critical_cols:
        n_missing = df[col].isna().sum()
        status = "✅" if n_missing == 0 else "❌"
        print(f"  {status} {col}: {n_missing} missing ({100*n_missing/len(df):.1f}%)")
        if n_missing > 0:
            no_missing = False
    
    if not no_missing:
        print("\n❌ FAILED: Missing values in critical columns!")
        return False
    
    # Check value ranges
    print("\n" + "=" * 80)
    print("VALUE RANGE CHECK")
    print("=" * 80)
    
    # Impact scores
    impact_min = df['impact_score'].min()
    impact_max = df['impact_score'].max()
    impact_valid = (impact_min >= 0) and (impact_max <= 1)
    status = "✅" if impact_valid else "❌"
    print(f"  {status} Impact scores: [{impact_min:.3f}, {impact_max:.3f}] (expected [0, 1])")
    
    # Allele frequencies
    freq_min = df['allele_frequency'].min()
    freq_max = df['allele_frequency'].max()
    freq_valid = (freq_min >= 0) and (freq_max <= 1)
    status = "✅" if freq_valid else "❌"
    print(f"  {status} Allele frequencies: [{freq_min:.6f}, {freq_max:.6f}] (expected [0, 1])")
    
    if not (impact_valid and freq_valid):
        print("\n❌ FAILED: Values out of valid range!")
        return False
    
    # Check distributions
    print("\n" + "=" * 80)
    print("DISTRIBUTION CHECK")
    print("=" * 80)
    
    # Impact source distribution
    print("\n  Impact Score Sources:")
    for source, count in df['impact_source'].value_counts().items():
        print(f"    {source:15s}: {count:5d} ({100*count/len(df):5.1f}%)")
    
    # Frequency source distribution
    print("\n  Allele Frequency Sources:")
    for source, count in df['frequency_source'].value_counts().items():
        print(f"    {source:15s}: {count:5d} ({100*count/len(df):5.1f}%)")
    
    # Protein distribution
    print("\n  Protein Distribution:")
    for protein, count in df['protein'].value_counts().items():
        print(f"    {protein:15s}: {count:5d} ({100*count/len(df):5.1f}%)")
    
    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    
    print("\n  Impact Scores:")
    print(f"    Mean:   {df['impact_score'].mean():.3f}")
    print(f"    Median: {df['impact_score'].median():.3f}")
    print(f"    Std:    {df['impact_score'].std():.3f}")
    
    print("\n  Allele Frequencies:")
    print(f"    Mean:   {df['allele_frequency'].mean():.6e}")
    print(f"    Median: {df['allele_frequency'].median():.6e}")
    print(f"    Std:    {df['allele_frequency'].std():.6e}")
    
    # Final checks
    print("\n" + "=" * 80)
    print("FINAL VALIDATION")
    print("=" * 80)
    
    checks = [
        ("All required columns present", all_present),
        ("No missing values in critical columns", no_missing),
        ("Impact scores in valid range [0, 1]", impact_valid),
        ("Allele frequencies in valid range [0, 1]", freq_valid),
        ("Total variants = 2139", len(df) == 2139),
    ]
    
    all_passed = True
    for check_name, passed in checks:
        status = "✅" if passed else "❌"
        print(f"  {status} {check_name}")
        if not passed:
            all_passed = False
    
    print("\n" + "=" * 80)
    if all_passed:
        print("✅ ✅ ✅  ALL CHECKS PASSED  ✅ ✅ ✅")
        print("=" * 80)
        print("\n🎉 Dataset is READY for Step 2.3: Domain → Trait-Bin Weight Assignment!")
        print("\nNext steps:")
        print("  1. Define trait bins (e.g., muscle strength, endurance, recovery)")
        print("  2. Assign domain-trait weights based on biological knowledge")
        print("  3. Create domain-trait mapping matrix")
        return True
    else:
        print("❌ ❌ ❌  VALIDATION FAILED  ❌ ❌ ❌")
        print("=" * 80)
        print("\n⚠️  Please fix the issues above before proceeding to Step 2.3")
        return False


if __name__ == '__main__':
    csv_file = Path('../protein_csvs/domain_variants_complete.csv')

    if not csv_file.exists():
        print(f"❌ Error: {csv_file} not found!")
        exit(1)

    success = verify_complete_dataset(csv_file)
    exit(0 if success else 1)

