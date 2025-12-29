# MechanobioAI: Genotype-to-Phenotype Prediction via Mechanobiology Knowledge Graph

**A mechanobiology-grounded synthetic model linking structural variation in titin-centered proteins to athletic performance traits.**

---

## Repository Structure

```
.
├── README.md                    # This file
├── domains_list.txt             # Simple list of all domains
│
├── proteins/                    # Raw protein data from UniProt
│   ├── {protein}.json           # UniProt JSON data
│   ├── {protein}_features.json  # Extracted features
│   ├── {protein}_domain_coordinates.txt
│   └── {protein}_variants.csv   # Per-protein variants
│
├── protein_csvs/                # Processed variant datasets
│   ├── all_variants.csv         # All variants (3,432)
│   ├── domain_variants.csv      # Domain variants only (2,139)
│   ├── domain_variants_with_impact.csv      # + impact scores
│   ├── domain_variants_complete.csv         # + allele frequencies (FINAL)
│   ├── IMPACT_SCORE_SUMMARY.md
│   └── COMPLETE_DATASET_SUMMARY.md
│
├── scripts/                     # All processing scripts
│   ├── build_variant_csv.py     # Step 1: Extract variants
│   ├── filter_domain_variants.py            # Step 2.1: Map to domains
│   ├── add_impact_scores.py                 # Step 2.2: Calculate impact
│   ├── add_allele_frequencies.py            # Step 2.2: Sample frequencies
│   ├── analyze_impact_scores.py             # Analysis
│   ├── analyze_allele_frequencies.py
│   ├── visualize_complete_dataset.py        # Visualizations
│   ├── visualize_final_impact_scores.py
│   └── verify_complete_dataset.py           # Validation
│
└── plots/                       # All visualizations
    ├── allele_frequency_distribution.png
    ├── complete_dataset_analysis.png
    ├── final_impact_score_analysis.png
    └── impact_score_distribution.png
```

**Note:** All scripts in `scripts/` should be run from within the `scripts/` directory. They use relative paths (`../protein_csvs/`, `../plots/`, `../proteins/`) to access data.

---

## Overview

This project builds a concrete, mechanobiology-grounded architecture demonstrating how observable athletic traits can carry recoverable signal about genotype-level variation in titin-related mechanosensing pathways.

### Core Concept

We link structural variation in a small titin-centered protein panel (**TTN, NBR1, p62/SQSTM1, MuRF1/TRIM63, MuRF2/TRIM55**) to three genomic burden traits (**power, endurance, recovery**) via a knowledge graph over:

```
Variants → Domains → Trait Bins → Phenotypes
```

From these burden traits, we generate synthetic observable phenotypes (strength, endurance, recovery capacity), then train a small ML model to invert this pipeline: **phenotype → genomic burden traits**.

---

## Architecture

### Conceptual Graph

```
       [Variant nodes]
      (per SNP in TTN, NBR1, p62, MuRF1, MuRF2)
                |
                | (in_domain)
                v
        [Domain nodes]
      (e.g. TTN_Kinase, NBR1_ZZ, p62_UBA, MuRF1_RING...)
                |
                | (relevance weights per bin)
                v
      [Genomic Trait Bin nodes]
   power_burden, endurance_burden, recovery_burden
                |
                | (linear generative model)
                v
        [Phenotype nodes]
   1RM, time_to_exhaustion, volume_tolerance, etc.
                ^
                |
      [Learned ML model fθ]
   (regressor/classifier: X_pheno → genomic trait bins)
```

**Key Components:**
- **Knowledge Graph (KG):** Variant → Domain → TraitBin
- **Generative Model:** TraitBins → Phenotypes
- **Trained Predictor:** Phenotypes → TraitBins (inverse mapping)

---

## Pipeline

### 1. Get Protein Mutation Data ✅ **COMPLETE**

**Goal:** Assemble a mechanobiology-relevant variant panel for our 5 proteins.

**Proteins:**
- **TTN** (Titin) - Sarcomeric spring, mechanosensor
  - *Note: For this MVP, only the Kinase domain is included. TTN is the largest known protein (~35,000 aa) with 100+ Ig/FN3 domains that would dominate the dataset. The Kinase domain is the mechanosensing signaling hub.*
- **NBR1** - Autophagy receptor
- **p62/SQSTM1** - Autophagy/ubiquitin adapter
- **MuRF1/TRIM63** - E3 ubiquitin ligase, muscle atrophy
- **MuRF2/TRIM55** - E3 ubiquitin ligase, muscle quality control

**Implementation:**
- Downloaded UniProt JSON for each protein
- Extracted all `VARIANT` features (position, alternative sequence, SIFT/PolyPhen predictions, population frequency)
- Extracted all `DOMAIN` features (begin, end, description)

**Output Files:**
- `proteins/{protein}.json` - Raw UniProt data
- `proteins/{protein}_features.json` - Extracted features
- `proteins/{protein}_domain_coordinates.txt` - Domain boundaries
- `proteins/{protein}_variants.csv` - Variant data per protein
- `all_variants.csv` - Combined variant table (3,432 variants)

**Variant Table Schema:**
```
protein, domain, position, wild_type, mutated_type, consequence_type,
sift_score, polyphen_score, allele_frequency
```

---

### 2. Build the Knowledge Graph: Variant → Domain → Trait Bins

#### 2.1 Map Variants to Domains ✅ **COMPLETE**

**Goal:** Assign each variant to its containing domain(s).

**Method:**
- For each variant at position `pos`, assign to all domains where `begin ≤ pos ≤ end`
- Variants outside annotated domains are excluded from the MVP

**Implementation:**
- Script: `filter_domain_variants.py`
- Filtered from 3,432 total variants → **2,139 domain variants**

**Output:**
- `domain_variants.csv` - Variants mapped to domains

**Domain Distribution:**
| Domain | Proteins | Variant Count |
|--------|----------|---------------|
| PB1 | p62, NBR1 | 419 |
| RING | MuRF1, MuRF2 | 318 |
| Ig | TTN | 289 |
| B-box | MuRF1, MuRF2 | 286 |
| COS | MuRF1, MuRF2 | 286 |
| Coiled-coil | MuRF1, MuRF2 | 286 |
| UBA | p62 | 215 |
| ZZ | NBR1 | 182 |
| FN3 | TTN | 120 |
| Kinase | TTN | 120 |

**Total: 10 unique domains**

---

#### 2.2 Turn Predictions/Frequencies into Numeric Values ✅ **COMPLETE**

**Goal:** Convert categorical predictions and missing values into numeric `impact_score` and `allele_frequency`.

##### Impact Score Calculation

**Hierarchical approach based on available data:**

| Source | Count | % | Formula |
|--------|-------|---|---------|
| Both SIFT & PolyPhen | 1,021 | 47.7% | `(1 - SIFT + PolyPhen) / 2` |
| SIFT only | 474 | 22.2% | `1 - SIFT` |
| PolyPhen only | 45 | 2.1% | `PolyPhen` |
| **Sampled** | **599** | **28.0%** | **Beta distribution by consequence type** |

**Sampling Strategy (for missing scores):**

| Consequence Type | Distribution | Mean | Rationale |
|-----------------|--------------|------|-----------|
| Stop-gain, Frameshift | Beta(8, 2) | ~0.80 | High impact (protein truncation) |
| Inframe deletion | Beta(6, 3) | ~0.67 | Moderate-high impact |
| Missense | Beta(3, 3) | ~0.50 | Variable impact |

**Impact Score Statistics:**
- Mean: 0.666
- Median: 0.670
- Range: [0.000, 1.000]
- Distribution: Bimodal (peaks at 0 and 1)

##### Allele Frequency Sampling

**Current Status:**
- Observed (from gnomAD): 1,145 (53.5%)
- Sampled: 994 (46.5%)

**Sampling Method:**
- **Distribution:** Log-Normal
- **Parameters:** μ = -5.217, σ = 0.924 (on log₁₀ scale)
- **Formula:** `10^(N(-5.217, 0.924))`
- **Rationale:** Best fit for rare variant frequencies, captures heavy tail

**Frequency Statistics:**
- Mean: 4.75 × 10⁻⁴
- Median: 4.76 × 10⁻⁶
- Range: [0.00, 0.412]

**Implementation:**
- Scripts: `add_impact_scores.py`, `add_allele_frequencies.py`
- Analysis: `analyze_impact_scores.py`, `analyze_allele_frequencies.py`
- Visualizations: `visualize_final_impact_scores.py`, `visualize_complete_dataset.py`

**Output:**
- `domain_variants_with_impact.csv` - Added impact scores
- `domain_variants_complete.csv` - **Final dataset with impact scores & frequencies**

**Validation:**
- ✅ Zero missing values in critical columns
- ✅ All scores in valid range [0, 1]
- ✅ Distributions match biological expectations
- ✅ Fully reproducible (seed = 42)

**Documentation:**
- `IMPACT_SCORE_SUMMARY.md` - Impact score methodology
- `COMPLETE_DATASET_SUMMARY.md` - Complete dataset documentation

---

#### 2.3 Domain → Trait-Bin Weights ⏳ **IN PROGRESS**

**Goal:** Define relevance weights for each domain to three genomic burden traits.

**Trait Bins:**
1. **power_burden** - Structural risk in power-relevant domains
2. **endurance_burden** - Endurance/fatigue-related pathways
3. **recovery_burden** - Repair/autophagy/ubiquitin pathways

**Method:**
For each domain, assign relevance scores [0, 2] for each trait bin based on:
- Protein function (from UniProt)
- Domain description
- Mechanobiology knowledge

**Target Output:**
```
domain_id, protein, power_weight, endurance_weight, recovery_weight
TTN_Kinase, TTN, 1.5, 0.5, 2.0
NBR1_ZZ, NBR1, 0.5, 0.5, 2.0
p62_UBA, p62, 0.0, 0.5, 2.0
MuRF1_RING, MuRF1, 1.0, 0.0, 2.0
...
```

**Status:** Ready to begin - see `domains_list.txt` for complete domain list

---

### 3. Mapping Genotypes → Genomic Trait Scores (via KG) ⏳ **PENDING**

**Goal:** Define how a genotype produces the 3 burden scores.

**For each individual:**

1. **Genotype dosage:** `g_i ∈ {0, 1, 2}` for each variant `i`

2. **Compute domain burdens:**
   ```
   domain_burden_d = Σ(impact_score_i × g_i) for all variants i in domain d
   ```

3. **Compute trait-bin scores:**
   ```
   y_t = Σ(domain_burden_d × w_d,t) for all domains d
   ```
   Where `w_d,t` comes from the domain weight table (Step 2.3)

**Output per individual:**
```
y = [y_power_burden, y_endurance_burden, y_recovery_burden]
```

---

### 4. Synthetic Genome Generation ⏳ **PENDING**

**Goal:** Create a synthetic cohort with realistic variation patterns.

**Method:**
- For each individual `n` and variant `i`:
  ```
  g_n,i ~ Binomial(2, p_i)
  ```
  where `p_i = allele_freq` from variant table

**Output:**
- **Genotype matrix:** `G ∈ ℝ^(N×V)` with entries in {0, 1, 2}
- **Domain burdens:** `B ∈ ℝ^(N×D)`
- **Trait-bin scores:** `Y ∈ ℝ^(N×3)` (power, endurance, recovery)

---

### 5. Mapping Genomic Trait Scores → Observable Phenotypes ⏳ **PENDING**

**Goal:** Generate plausible phenotypes from the 3 burden traits.

**Phenotype Vector (per individual):**
- `strength_1RM` (normalized)
- `time_to_exhaustion`
- `volume_tolerance` (max sustainable weekly volume)
- `muscle_mass_index`
- `resting_HR` (or HRR)

**Generative Model:**
```
x_pheno,n = A · y_n + b + ε_n
```

Where:
- `y_n = [y_power_burden, y_endurance_burden, y_recovery_burden]ᵀ`
- `A` is a hand-chosen matrix (more burden = worse performance):
  ```
  # rows: phenotypes; cols: [power, endurance, recovery]
  A = [
    [-1.2, -0.3, -0.4],   # strength_1RM
    [-0.2, -1.5, -0.5],   # time_to_exhaustion
    [-0.3, -0.4, -1.5],   # volume_tolerance
    [-1.0, -0.6, -0.3],   # muscle_mass_index
    [ 0.0,  1.2,  0.5],   # resting_HR (higher burden → higher HR)
  ]
  ```
- `b` is baseline mean vector
- `ε_n` is Gaussian noise

**Output:**
- **Phenotype matrix:** `X_pheno ∈ ℝ^(N×P)` (P = number of phenotypes)
- **Dataset:** Pairs `(X_pheno, Y)`

---

### 6. Train Prediction Model: Phenotypes → Genomic Trait Bins ⏳ **PENDING**

**Goal:** Demonstrate that observable traits carry recoverable signal about genomic burden.

#### 6.1 Regression Model

**Inputs:** `X_pheno` (e.g., 5 phenotypes)
**Targets:** `Y` (3 continuous burden scores)

**Model Options:**
- Linear regressor
- Small MLP with 1-2 hidden layers (16-32 units)

**Loss:** MSE on 3 outputs

**Evaluation:**
- R² per trait
- Predicted vs true plots for each burden trait

#### 6.2 Classification Model (Optional)

**Method:**
- Discretize each `y_t`: "high burden" (top 30%), "low burden" (bottom 30%)
- Train binary classifiers from `X_pheno` to labels

**Loss:** Cross-entropy
**Metrics:** Accuracy, AUROC per trait

---



## Current Status

### ✅ Completed

- [x] **Step 1:** Protein mutation data extraction
  - 5 proteins, 3,432 total variants
  - 10 unique domains identified

- [x] **Step 2.1:** Variant → Domain mapping
  - 2,139 domain variants

- [x] **Step 2.2:** Impact scores & allele frequencies
  - 100% complete, zero missing values
  - Biologically realistic distributions
  - Fully reproducible (seed = 42)

### ⏳ In Progress

- [ ] **Step 2.3:** Domain → Trait-Bin weights
  - Define power/endurance/recovery weights for 10 domains

### 📋 Pending

- [ ] **Step 3:** Genotype → Genomic trait score mapping
- [ ] **Step 4:** Synthetic genome generation
- [ ] **Step 5:** Genomic traits → Phenotype generation
- [ ] **Step 6:** Train inverse model (Phenotype → Genomic traits)

---

## Key Results (So Far)

### Dataset Statistics

| Metric | Value |
|--------|-------|
| Total variants | 2,139 |
| Unique domains | 10 |
| Proteins | 5 |
| Impact score mean | 0.666 |
| Allele frequency median | 4.76 × 10⁻⁶ |
| Missing values | 0 |

### Protein Distribution

| Protein | Variants | % |
|---------|----------|---|
| p62 | 634 | 29.6% |
| MuRF2 | 604 | 28.2% |
| TTN | 409 | 19.1% |
| MuRF1 | 310 | 14.5% |
| NBR1 | 182 | 8.5% |

### Visualizations

- `allele_frequency_distribution.png` - Frequency analysis (6 panels)
- `complete_dataset_analysis.png` - Comprehensive overview (8 panels)
- `final_impact_score_analysis.png` - Impact score distributions

---

## Next Steps

1. **Define domain-trait weights** (Step 2.3)
   - Use `domains_list.txt` as reference
   - Assign power/endurance/recovery weights [0, 2] for each domain
   - Create `domain_trait_weights.csv`

2. **Implement genotype-to-burden mapping** (Step 3)

3. **Generate synthetic cohort** (Step 4)

4. **Create phenotype generation model** (Step 5)

5. **Train and evaluate inverse model** (Step 6)

---

## Reproducibility

All analyses are fully reproducible:
- **Random seed:** 42
- **Python version:** 3.x
- **Key packages:** pandas, numpy, scipy, matplotlib, seaborn

**To reproduce Steps 1-2.2:**
```bash
cd scripts
python build_variant_csv.py
python filter_domain_variants.py
python add_impact_scores.py
python add_allele_frequencies.py
python verify_complete_dataset.py
```

**Note:** All scripts must be run from the `scripts/` directory as they use relative paths to access data files.

---

## References

- **UniProt:** Protein data source
- **gnomAD:** Population allele frequencies
- **SIFT/PolyPhen:** Variant effect predictions

---

## License

Research project - AAAI-IAAA 2025

---

**Last Updated:** 2025-12-29
**Status:** Step 2.2 Complete, Ready for Step 2.3

