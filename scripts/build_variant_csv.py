#!/usr/bin/env python3
"""
Build a CSV file with all variants from protein JSON files.
Maps variants to domains and extracts prediction scores and allele frequencies.
"""

import json
import csv
import sys
from pathlib import Path
from typing import List, Dict, Any, Tuple
import numpy as np


def parse_domain_coordinates(domain_file: Path) -> List[Dict[str, Any]]:
    """Parse domain coordinates from text file."""
    domains = []
    with open(domain_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Parse format: "ZZ = [191, 243] (Domain, ZZ-type)"
            parts = line.split('=')
            if len(parts) != 2:
                continue
            
            domain_name = parts[0].strip()
            rest = parts[1].strip()
            
            # Extract coordinates
            coord_start = rest.find('[')
            coord_end = rest.find(']')
            if coord_start == -1 or coord_end == -1:
                continue
            
            coords = rest[coord_start+1:coord_end].split(',')
            begin = int(coords[0].strip())
            end = int(coords[1].strip())
            
            # Extract description (everything in parentheses)
            desc_start = rest.find('(')
            desc_end = rest.find(')')
            description = rest[desc_start+1:desc_end].strip() if desc_start != -1 else ""
            
            domains.append({
                'name': domain_name,
                'begin': begin,
                'end': end,
                'description': description
            })
    
    return domains


def find_domains_for_position(position: int, domains: List[Dict[str, Any]]) -> List[str]:
    """Find all domains that contain the given position."""
    matching_domains = []
    for domain in domains:
        if domain['begin'] <= position <= domain['end']:
            matching_domains.append(domain['name'])
    return matching_domains


def extract_variant_data(variant: Dict[str, Any], protein: str, domains: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Extract variant data and map to domains."""
    position = int(variant.get('begin', 0))
    
    # Find matching domains
    matching_domains = find_domains_for_position(position, domains)
    
    # If no domain match, use "None" as domain
    if not matching_domains:
        matching_domains = ["None"]
    
    # Extract SIFT and PolyPhen scores
    sift_score = np.nan
    polyphen_score = np.nan
    
    if 'predictions' in variant:
        for pred in variant['predictions']:
            algo = pred.get('predAlgorithmNameType', '').upper()
            score = pred.get('score')
            if algo == 'SIFT' and score is not None:
                sift_score = score
            elif algo == 'POLYPHEN' and score is not None:
                polyphen_score = score
    
    # Extract allele frequency (first available)
    allele_frequency = np.nan
    if 'populationFrequencies' in variant:
        pop_freqs = variant['populationFrequencies']
        if pop_freqs and len(pop_freqs) > 0:
            allele_frequency = pop_freqs[0].get('frequency', np.nan)
    
    # Extract other fields
    wild_type = variant.get('wildType', '')
    mutated_type = variant.get('mutatedType', '')
    consequence_type = variant.get('consequenceType', '')
    
    # Create one row per domain
    rows = []
    for domain in matching_domains:
        rows.append({
            'protein': protein,
            'domain': domain,
            'position': position,
            'wild_type': wild_type,
            'mutated_type': mutated_type,
            'consequence_type': consequence_type,
            'sift_score': sift_score,
            'polyphen_score': polyphen_score,
            'allele_frequency': allele_frequency
        })
    
    return rows


def build_variant_csv(protein: str, json_file: Path, domain_file: Path, output_file: Path):
    """Build CSV file with variants mapped to domains."""

    # Load domain coordinates
    print(f"Loading domain coordinates from {domain_file}...")
    domains = parse_domain_coordinates(domain_file)
    print(f"  Found {len(domains)} domains: {[d['name'] for d in domains]}")

    # Load variant JSON
    print(f"Loading variants from {json_file}...")
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Extract all variants
    features = data.get('features', [])
    variants = [f for f in features if f.get('type') == 'VARIANT']
    print(f"  Found {len(variants)} variants")

    # Process variants
    all_rows = []
    multi_domain_count = 0
    no_domain_count = 0
    domain_with_scores_count = 0

    for variant in variants:
        rows = extract_variant_data(variant, protein, domains)

        if len(rows) > 1:
            multi_domain_count += 1
        if rows and rows[0]['domain'] == "None":
            no_domain_count += 1

        # Check if variant is in a domain AND has prediction scores
        for row in rows:
            if row['domain'] != "None":
                has_sift = not (isinstance(row['sift_score'], float) and np.isnan(row['sift_score']))
                has_polyphen = not (isinstance(row['polyphen_score'], float) and np.isnan(row['polyphen_score']))
                if has_sift or has_polyphen:
                    domain_with_scores_count += 1
                    break  # Count each variant only once

        all_rows.extend(rows)

    # Write CSV
    print(f"Writing {len(all_rows)} rows to {output_file}...")
    fieldnames = ['protein', 'domain', 'position', 'wild_type', 'mutated_type',
                  'consequence_type', 'sift_score', 'polyphen_score', 'allele_frequency']

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    # Summary
    print(f"\nSummary:")
    print(f"  Total variants: {len(variants)}")
    print(f"  Total rows (variant-domain pairs): {len(all_rows)}")
    print(f"  Variants in multiple domains: {multi_domain_count}")
    print(f"  Variants with no domain match: {no_domain_count}")
    print(f"  Variants in domain WITH SIFT or PolyPhen score: {domain_with_scores_count}")
    print(f"\nDone!")


if __name__ == '__main__':
    # Define all proteins to process
    proteins_dir = Path('proteins')
    proteins = [
        {
            'name': 'NBR1',
            'json': proteins_dir / 'nbr1.json',
            'domains': proteins_dir / 'nbr1_domain_coordinates.txt'
        },
        {
            'name': 'p62',
            'json': proteins_dir / 'p62.json',
            'domains': proteins_dir / 'p62_domain_coordinates.txt'
        },
        {
            'name': 'TTN',
            'json': proteins_dir / 'titin.json',
            'domains': proteins_dir / 'titin_domain_coordinates.txt'
        },
        {
            'name': 'MuRF1',
            'json': proteins_dir / 'murf1.json',
            'domains': proteins_dir / 'murf1_domain_coordinates.txt'
        },
        {
            'name': 'MuRF2',
            'json': proteins_dir / 'murf2.json',
            'domains': proteins_dir / 'murf2_domain_coordinates.txt'
        }
    ]

    print("=" * 80)
    print("Building combined variant CSV file for all proteins")
    print("=" * 80)

    all_rows = []
    total_stats = {
        'total_variants': 0,
        'multi_domain_count': 0,
        'no_domain_count': 0,
        'domain_with_scores_count': 0
    }

    for protein_info in proteins:
        print(f"\n{'=' * 80}")
        print(f"Processing {protein_info['name']}")
        print(f"{'=' * 80}")

        json_file = protein_info['json']
        domain_file = protein_info['domains']

        # Check if files exist
        if not json_file.exists():
            print(f"  ⚠️  WARNING: {json_file} not found, skipping...")
            continue
        if not domain_file.exists():
            print(f"  ⚠️  WARNING: {domain_file} not found, skipping...")
            continue

        # Load domain coordinates
        print(f"Loading domain coordinates from {domain_file}...")
        domains = parse_domain_coordinates(domain_file)
        print(f"  Found {len(domains)} domains: {[d['name'] for d in domains]}")

        # Load variant JSON
        print(f"Loading variants from {json_file}...")
        with open(json_file, 'r') as f:
            data = json.load(f)

        # Extract all variants
        features = data.get('features', [])
        variants = [f for f in features if f.get('type') == 'VARIANT']
        print(f"  Found {len(variants)} variants")

        # Process variants
        protein_rows = []
        multi_domain_count = 0
        no_domain_count = 0
        domain_with_scores_count = 0

        for variant in variants:
            rows = extract_variant_data(variant, protein_info['name'], domains)

            if len(rows) > 1:
                multi_domain_count += 1
            if rows and rows[0]['domain'] == "None":
                no_domain_count += 1

            # Check if variant is in a domain AND has prediction scores
            for row in rows:
                if row['domain'] != "None":
                    has_sift = not (isinstance(row['sift_score'], float) and np.isnan(row['sift_score']))
                    has_polyphen = not (isinstance(row['polyphen_score'], float) and np.isnan(row['polyphen_score']))
                    if has_sift or has_polyphen:
                        domain_with_scores_count += 1
                        break  # Count each variant only once

            protein_rows.extend(rows)

        # Add to combined list
        all_rows.extend(protein_rows)

        # Update total stats
        total_stats['total_variants'] += len(variants)
        total_stats['multi_domain_count'] += multi_domain_count
        total_stats['no_domain_count'] += no_domain_count
        total_stats['domain_with_scores_count'] += domain_with_scores_count

        # Summary for this protein
        print(f"\nSummary for {protein_info['name']}:")
        print(f"  Total variants: {len(variants)}")
        print(f"  Total rows (variant-domain pairs): {len(protein_rows)}")
        print(f"  Variants in multiple domains: {multi_domain_count}")
        print(f"  Variants with no domain match: {no_domain_count}")
        print(f"  Variants in domain WITH SIFT or PolyPhen score: {domain_with_scores_count}")

    # Write combined CSV
    output_file = Path('../protein_csvs/all_variants.csv')
    print(f"\n{'=' * 80}")
    print(f"Writing combined CSV with {len(all_rows)} rows to {output_file}...")
    print(f"{'=' * 80}")

    fieldnames = ['protein', 'domain', 'position', 'wild_type', 'mutated_type',
                  'consequence_type', 'sift_score', 'polyphen_score', 'allele_frequency']

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    # Overall summary
    print(f"\n{'=' * 80}")
    print("OVERALL SUMMARY")
    print(f"{'=' * 80}")
    print(f"  Total variants across all proteins: {total_stats['total_variants']}")
    print(f"  Total rows in combined CSV: {len(all_rows)}")
    print(f"  Variants in multiple domains: {total_stats['multi_domain_count']}")
    print(f"  Variants with no domain match: {total_stats['no_domain_count']}")
    print(f"  Variants in domain WITH SIFT or PolyPhen score: {total_stats['domain_with_scores_count']}")
    print(f"\n✅ Combined CSV saved to: {output_file}")
    print(f"{'=' * 80}")

