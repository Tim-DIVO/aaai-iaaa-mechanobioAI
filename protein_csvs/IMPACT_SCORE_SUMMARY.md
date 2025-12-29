# Impact Score Implementation Summary

## Overview
Successfully implemented step 2.2 of the MVP: converting SIFT/PolyPhen predictions into unified impact scores for all domain variants.

## Files Created

### Data Files
- **`domain_variants.csv`** - Filtered variants (only those in protein domains): 2,139 variants
- **`domain_variants_with_impact.csv`** - Final dataset with impact scores: 2,139 variants

### Scripts
- **`filter_domain_variants.py`** - Filters variants to only those in domains
- **`add_impact_scores.py`** - Adds impact scores to variant data
- **`analyze_impact_scores.py`** - Analyzes impact score distributions
- **`visualize_final_impact_scores.py`** - Creates comprehensive visualizations

### Visualizations
- **`impact_score_distribution.png`** - Initial analysis of available scores
- **`final_impact_score_analysis.png`** - Comprehensive final analysis

## Impact Score Calculation Method

### 1. When Both SIFT and PolyPhen Available (47.7% of variants)
```
impact_from_sift = 1 - sift_score
impact_from_polyphen = polyphen_score
impact_score = (impact_from_sift + impact_from_polyphen) / 2
```

### 2. When Only SIFT Available (22.2% of variants)
```
impact_score = 1 - sift_score
```

### 3. When Only PolyPhen Available (2.1% of variants)
```
impact_score = polyphen_score
```

### 4. When Neither Available (28.0% of variants)
Sample from Beta distribution based on consequence type:

| Consequence Type | Beta(a, b) | Mean | Rationale |
|-----------------|------------|------|-----------|
| `stop gained` | Beta(8, 2) | ~0.80 | Loss of function, highly damaging |
| `frameshift` | Beta(8, 2) | ~0.80 | Disrupts reading frame, highly damaging |
| `initiator codon variant` | Beta(8, 2) | ~0.80 | Affects start codon, very damaging |
| `inframe deletion` | Beta(6, 3) | ~0.67 | Removes amino acids, likely damaging |
| `delins` | Beta(6, 3) | ~0.67 | Deletion+insertion, likely damaging |
| `missense` | Beta(3, 3) | ~0.50 | Variable impact, centered distribution |
| `insertion` | Beta(3, 3) | ~0.50 | Moderate impact |
| Unknown (`-`) | Beta(2, 2) | ~0.50 | No information |

## Results Summary

### Overall Statistics (n=2,139)
- **Mean**: 0.666
- **Median**: 0.670
- **Std Dev**: 0.271
- **Range**: [0.000, 1.000]

### By Source
| Source | Count | % | Mean | Median | Std |
|--------|-------|---|------|--------|-----|
| Both scores | 1,021 | 47.7% | 0.628 | 0.566 | 0.269 |
| SIFT only | 474 | 22.2% | 0.860 | 0.960 | 0.231 |
| PolyPhen only | 45 | 2.1% | 0.690 | 0.938 | 0.377 |
| Sampled | 599 | 28.0% | 0.575 | 0.579 | 0.213 |

### By Protein
| Protein | Variants | Mean Impact Score |
|---------|----------|-------------------|
| p62 | 634 | ~0.70 |
| MuRF2 | 604 | ~0.68 |
| TTN | 409 | ~0.65 |
| MuRF1 | 310 | ~0.66 |
| NBR1 | 182 | ~0.64 |

### Sampled Scores by Consequence Type
| Consequence Type | Count | Mean Impact |
|-----------------|-------|-------------|
| `stop gained` | 86 | 0.822 |
| `frameshift` | 49 | 0.798 |
| `initiator codon variant` | 1 | 0.857 |
| `inframe deletion` | 14 | 0.677 |
| `missense` | 440 | 0.500 |
| `insertion` | 3 | 0.542 |
| `delins` | 1 | 0.509 |
| `-` | 5 | 0.439 |

## Key Observations

1. **Bimodal Distribution**: Impact scores show peaks near 0 (benign) and 1 (damaging), reflecting the biological reality of variant effects.

2. **SIFT Bias**: SIFT-only variants are heavily skewed toward high impact (mean 0.86), suggesting SIFT may be more conservative in its predictions.

3. **Balanced Overall**: The combination of all sources yields a balanced distribution (mean 0.666, median 0.670).

4. **Sampled Scores Match Expectations**: Sampled scores for different consequence types align well with biological expectations:
   - Stop-gain/frameshift: ~0.80 (high impact)
   - Inframe deletions: ~0.67 (moderate-high)
   - Missense: ~0.50 (variable)

5. **Reproducibility**: All sampling uses seed=42 for reproducible results.

## Next Steps (MVP Step 2.3)

The next step is to create domain → trait-bin weights:
- Define relevance scores for each domain to power/endurance/recovery burden
- Use LLM to assign weights based on protein function and domain description
- Create domain weight table with columns: `domain_id, protein, power_weight, endurance_weight, recovery_weight`

## Data Schema

### Final CSV Columns
```
protein              : Protein name (NBR1, p62, TTN, MuRF1, MuRF2)
domain               : Domain name (ZZ, PB1, UBA, RING, Kinase, etc.)
position             : Amino acid position
wild_type            : Original amino acid
mutated_type         : Mutated amino acid
consequence_type     : Type of mutation (missense, stop gained, etc.)
sift_score           : Original SIFT score (0-1, lower=damaging)
polyphen_score       : Original PolyPhen score (0-1, higher=damaging)
allele_frequency     : Population allele frequency
impact_score         : Unified impact score (0-1, higher=damaging)
impact_source        : Source of impact score (both/sift_only/polyphen_only/sampled)
```

## Validation

✅ All variants have impact scores (no missing values)  
✅ Impact scores are in valid range [0, 1]  
✅ Distribution matches biological expectations  
✅ Sampled scores align with consequence severity  
✅ Reproducible (seed=42)  
✅ Well-documented and visualized

