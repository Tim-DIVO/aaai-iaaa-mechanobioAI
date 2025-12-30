# MechanobioAI: Genotype-to-Phenotype Prediction via Mechanobiology Knowledge Graph

**A mechanobiology-grounded synthetic model linking structural variation in titin-centered proteins to athletic performance traits.**

---

## Repository Structure

```
.
├── README.md                    # This file
├── METHODS_RESULTS.md           # Paper-style methods, results, discussion
├── abstract.txt                 # Conference abstract
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
│   └── domain_variants_complete.csv         # + allele frequencies (FINAL)
│
├── trait_mapping/               # Domain → Trait weights
│   ├── combined_domain.csv      # Final weights (Opus + Chat averaged)
│   ├── phenotype_definitions.csv # 10 phenotype coefficients
│   ├── opus_domain.csv          # Claude Opus weights
│   └── chat_domain.csv          # ChatGPT-4o weights
│
├── synthetic_cohort/            # Generated synthetic population
│   ├── synthetic_cohort_summary.csv  # 50k individuals: burdens + 10 phenotypes
│   ├── synthetic_genotypes.csv       # Sparse genotype matrix
│   ├── model_results.json            # Model performance metrics
│   └── *.npy                         # Saved predictions for plotting
│
├── scripts/                     # All processing scripts
│   ├── build_variant_csv.py              # Step 1: Extract variants
│   ├── filter_domain_variants.py         # Step 2.1: Map to domains
│   ├── add_impact_scores.py              # Step 2.2: Calculate impact
│   ├── add_allele_frequencies.py         # Step 2.2: Sample frequencies
│   ├── generate_synthetic_cohort.py      # Step 3-4: Generate cohort
│   ├── train_burden_predictor.py         # Step 5: Train MLP model
│   ├── create_analysis_plots.py          # Variant & cohort plots
│   └── create_model_plots.py             # Model evaluation plots
│
└── plots/                       # All visualizations (17 figures)
    ├── variant_*.png            # Variant analysis plots
    ├── cohort_*.png             # Synthetic cohort plots
    ├── model_*.png              # Model evaluation plots
    └── fig*.png                 # Poster figures
```

**Note:** All scripts in `scripts/` should be run from within the `scripts/` directory. They use relative paths to access data.

---

## Overview

This project builds a concrete, mechanobiology-grounded architecture demonstrating how observable athletic traits can carry recoverable signal about genotype-level variation in titin-related mechanosensing pathways.

### Core Concept

We link structural variation in a small titin-centered protein panel (**TTN, NBR1, p62/SQSTM1, MuRF1/TRIM63, MuRF2/TRIM55**) to two pathway burden scores (**hypertrophy, autophagy**) via a knowledge graph over:

```
Variants → Domains → Pathway Burdens → 10 Phenotypes
```

From these burden scores, we generate 10 observable phenotypes (muscle mass, strength, VO2max, recovery, etc.), then train a neural network to invert this pipeline: **10 phenotypes → 2 pathway burdens** with R² > 0.99.

---

## Architecture

### Conceptual Graph

```
       [Variant nodes]
      (2,139 SNPs in TTN, NBR1, p62, MuRF1, MuRF2)
                |
                | (in_domain, impact_score)
                v
        [Domain nodes]
      (21 domains: TTN_Kinase, NBR1_ZZ, p62_UBA, MuRF1_RING...)
                |
                | (pathway weights from literature)
                v
      [Pathway Burden nodes]
   hypertrophy_burden, autophagy_burden
                |
                | (linear generative model with pathway-specific coefficients)
                v
        [Phenotype nodes]
   muscle_growth_rate, natural_muscle_mass, strength_1rm,
   power_output, vo2max, time_to_exhaustion, recovery_speed,
   doms_severity, injury_susceptibility, training_responsiveness
                ^
                |
      [Learned MLP model fθ]
   (regressor: 10 phenotypes → 2 pathway burdens, R² > 0.99)
```

**Key Components:**
- **Knowledge Graph (KG):** Variant → Domain → Pathway (with literature-derived weights)
- **Generative Model:** Pathway Burdens → 10 Phenotypes (with pathway-specific coefficients)
- **Trained Predictor:** 10 Phenotypes → 2 Pathway Burdens (inverse mapping, R² > 0.99)

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

**Domain-Protein Distribution (21 combinations):**
| Domain | Protein | Variant Count |
|--------|---------|---------------|
| Kinase | TTN | 409 |
| IDR | MuRF2 | 284 |
| IDR2 | p62 | 242 |
| PB1 | p62 | 158 |
| Coiled-coil | MuRF2 | 116 |
| COS | MuRF2 | 85 |
| Coiled-coil | MuRF1 | 81 |
| UBA | p62 | 77 |
| ZZ | p62 | 76 |
| RING | MuRF1 | 69 |
| IDR1 | p62 | 69 |
| COS | MuRF1 | 68 |
| RING | MuRF2 | 67 |
| B-box | MuRF1 | 66 |
| B-box | MuRF2 | 52 |
| ZZ | NBR1 | 50 |
| IDR2 | NBR1 | 48 |
| Coiled-coil | NBR1 | 47 |
| IDR1 | NBR1 | 37 |
| IDR | MuRF1 | 26 |
| LIR | p62 | 12 |

**Total: 21 domain-protein combinations** (see `domains_list.txt` for details)

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

#### 2.3 Domain → Trait Weights ✅ **COMPLETE**

**Goal:** Define relevance weights for each domain to two biological pathways.

**Trait Pathways:**
1. **Hypertrophy** - Anabolic signaling, muscle growth
2. **Autophagy** - Catabolic signaling, protein turnover

**Method:**
- LLM-assisted literature review (Opus + ChatGPT-4o, averaged)
- Weights assigned based on domain involvement in each pathway [0, 1]
- References from primary literature (Lange 2005, Bogomolovas 2021, etc.)

**Output:**
- `trait_mapping/combined_domain.csv` - 21 domain-protein combinations with hypertrophy/autophagy weights

**Key Insight - Pathway Directionality:**
- **Hypertrophy** = anabolic (builds muscle) → Disruption is BAD
- **Autophagy** = catabolic (breaks down proteins) → Disruption is GOOD (less catabolism)

---

### 3. Mapping Genotypes → Burden Scores → Phenotypes ✅ **COMPLETE**

**Goal:** Define how a genotype produces burden scores and 10 observable phenotypes.

**For each individual:**

1. **Genotype dosage:** `g_i ∈ {0, 1, 2}` for each variant `i`

2. **Compute pathway burdens:**
   ```
   hypertrophy_burden = Σ(impact_score_i × hypertrophy_weight_i × g_i)
   autophagy_burden = Σ(impact_score_i × autophagy_weight_i × g_i)
   ```

3. **Compute 10 observable phenotypes:**
   ```
   phenotype = baseline + hyp_coef × hyp_burden + auto_coef × auto_burden + noise
   ```

**Phenotype Definitions (from `trait_mapping/phenotype_definitions.csv`):**

| Phenotype | Hyp Coef | Auto Coef | Dominant Pathway |
|-----------|----------|-----------|------------------|
| muscle_growth_rate | -1.5 | +0.5 | Hypertrophy |
| natural_muscle_mass | -1.5 | +0.6 | Hypertrophy |
| strength_1rm | -1.2 | +0.3 | Hypertrophy |
| power_output | -1.3 | +0.2 | Hypertrophy |
| vo2max | -0.2 | -1.5 | Autophagy |
| time_to_exhaustion | -0.1 | -1.5 | Autophagy |
| recovery_speed | -0.2 | -1.5 | Autophagy |
| doms_severity | +0.2 | +1.5 | Autophagy |
| injury_susceptibility | +0.4 | +0.8 | Both |
| training_responsiveness | -0.6 | -0.6 | Both |

**Key Insight - Pathway Effects:**
- **Hypertrophy disruption** = can't build muscle = bad for strength/power phenotypes
- **Autophagy disruption** = less catabolism = GOOD for muscle mass, but BAD for endurance/recovery (impaired mitochondrial quality control)

---

### 4. Synthetic Genome Generation ✅ **COMPLETE**

**Goal:** Create a synthetic cohort with realistic variation patterns.

**Method:**
- For each individual `n` and variant `i`:
  ```
  g_n,i ~ Binomial(2, p_i)
  ```
  where `p_i = allele_freq` from variant table
- Script: `scripts/generate_synthetic_cohort.py`

**Output (in `synthetic_cohort/`):**
- `synthetic_cohort_summary.csv` - 50,000 individuals with burden scores + 10 phenotypes
- `synthetic_genotypes.csv` - Sparse genotype matrix (only non-zero entries)

**Statistics (N=50,000):**
| Metric | Mean | Std |
|--------|------|-----|
| Hypertrophy burden | 0.434 | 0.307 |
| Autophagy burden | 0.578 | 0.385 |
| Variants per individual | 1.78 | 1.05 |

**Partial Correlations (unique pathway contributions):**

| Phenotype | r(Hyp\|Auto) | r(Auto\|Hyp) | Dominant |
|-----------|--------------|--------------|----------|
| muscle_growth_rate | -0.87 | +0.60 | Hypertrophy |
| natural_muscle_mass | -0.91 | +0.73 | Hypertrophy |
| strength_1rm | -0.75 | +0.34 | Hypertrophy |
| power_output | -0.84 | +0.28 | Hypertrophy |
| vo2max | -0.19 | -0.87 | Autophagy |
| time_to_exhaustion | -0.12 | -0.91 | Autophagy |
| recovery_speed | -0.23 | -0.91 | Autophagy |
| doms_severity | +0.18 | +0.87 | Autophagy |
| training_responsiveness | -0.58 | -0.67 | Both |

---

### 5. Train Prediction Model: Phenotypes → Burden Scores ✅ **COMPLETE**

**Goal:** Demonstrate that 10 observable phenotypes carry recoverable signal about 2 latent genomic burdens.

**Regression Model:**
- **Input:** 10 phenotypes (muscle_growth_rate, natural_muscle_mass, ..., training_responsiveness)
- **Targets:** `[hypertrophy_burden, autophagy_burden]` (2 latent scores)

**Model Architecture:**
- MLP with 2 hidden layers (64, 32 neurons)
- ReLU activation, Adam optimizer
- Early stopping with 10% validation split
- 85/15 train/test split (42,500 / 7,500 samples)

**Results:**

| Target | Train R² | Test R² | Test MSE | Test MAE |
|--------|----------|---------|----------|----------|
| Hypertrophy Burden | 0.9910 | 0.9902 | 0.00088 | 0.0187 |
| Autophagy Burden | 0.9939 | 0.9935 | 0.00092 | 0.0195 |

**Key Finding:** Near-perfect recovery (R² > 0.99) of latent genomic burden scores from 10 observable phenotypes, demonstrating that phenotype patterns carry strong signal about underlying genetic pathway disruption.

---

## Current Status

### ✅ Completed

- [x] **Step 1:** Protein mutation data extraction
  - 5 proteins, 3,432 total variants
  - 21 domain-protein combinations identified

- [x] **Step 2.1:** Variant → Domain mapping
  - 2,139 domain variants

- [x] **Step 2.2:** Impact scores & allele frequencies
  - 100% complete, zero missing values
  - Biologically realistic distributions
  - Fully reproducible (seed = 42)

- [x] **Step 2.3:** Domain → Trait weights
  - Hypertrophy/autophagy weights for 21 domain-protein combinations
  - LLM-assisted literature review (Opus + ChatGPT-4o averaged)

- [x] **Step 3-4:** Synthetic cohort generation
  - 50,000 synthetic individuals
  - 10 observable phenotypes with pathway-specific coefficients
  - Genotype → Burden → Phenotypes pipeline implemented

- [x] **Step 5:** Inverse prediction model
  - MLP (64-32) trained with early stopping
  - R² > 0.99 on both burden scores
  - No overfitting (train ≈ test performance)

---

## Key Results

### Dataset Statistics

| Metric | Value |
|--------|-------|
| Total variants | 2,139 |
| Domain-protein combinations | 21 |
| Proteins | 5 |
| Observable phenotypes | 10 |
| Impact score mean | 0.666 |
| Allele frequency median | 4.76 × 10⁻⁶ |
| Synthetic individuals | 50,000 |

### Model Performance

| Burden Score | Test R² | Interpretation |
|--------------|---------|----------------|
| Hypertrophy | 0.990 | Near-perfect recovery |
| Autophagy | 0.994 | Near-perfect recovery |

### Protein Distribution

| Protein | Variants | % |
|---------|----------|---|
| p62 | 634 | 29.6% |
| MuRF2 | 604 | 28.2% |
| TTN | 409 | 19.1% |
| MuRF1 | 310 | 14.5% |
| NBR1 | 182 | 8.5% |

### Visualizations

**Variant Analysis:**
- `plots/variant_impact_scores.png` - Impact scores (observed vs imputed)
- `plots/variant_allele_frequencies.png` - Allele frequencies by source and protein

**Synthetic Cohort:**
- `plots/cohort_phenotype_distributions.png` - 10 phenotype distributions
- `plots/cohort_burden_phenotype_scatter.png` - Burden scores colored by phenotypes

**Model Evaluation:**
- `plots/model_training_loss.png` - Training loss curve
- `plots/model_pred_vs_true.png` - Predicted vs true burden scores
- `plots/model_error_distribution.png` - Prediction error histograms
- `plots/model_performance_summary.png` - Train vs test R² comparison

---

## Pipeline Complete ✅

All 5 steps of the genotype-to-phenotype pipeline have been implemented and validated:

1. ✅ Variant extraction from UniProt
2. ✅ Domain mapping and impact scoring
3. ✅ Pathway weight assignment (LLM-assisted)
4. ✅ Synthetic cohort generation (50k individuals, 10 phenotypes)
5. ✅ Inverse prediction model (R² > 0.99)

---

## Conceptual Limitations

Below are the things that would turn this into an actual working model:

1. **Biophysical modeling for each mutation** to figure out changes in binding strength / protein dynamics in relevant pathways.
   This requires:
   - a) Accurate structure prediction for proteins altered by single-point, non-deleterious mutations (an issue especially for huge proteins like Titin, so currently not technically possible)
   - b) Run a whole bunch of MD simulations / binding predictions (e.g., AlphaFold/AlphaProteo) to figure out how these new structures behave

2. **Extract actual knowledge graph ontologies** from all of PubMed to link proteins (and mutations therein) to pathways and phenotypes in an exhaustive manner.

3. **Get actual phenotype and genotype data**, ideally alongside training logs, to validate biophysical and KG models.

---

## Reproducibility

All analyses are fully reproducible:
- **Random seed:** 42
- **Python version:** 3.x
- **Key packages:** pandas, numpy, scipy, matplotlib, seaborn

**To reproduce full pipeline:**
```bash
cd scripts
python build_variant_csv.py
python filter_domain_variants.py
python add_impact_scores.py
python add_allele_frequencies.py
python generate_synthetic_cohort.py
python train_burden_predictor.py
python create_analysis_plots.py
python create_model_plots.py
```

**Note:** All scripts must be run from the `scripts/` directory as they use relative paths.

---

## References

- **UniProt:** Protein data source
- **gnomAD:** Population allele frequencies
- **SIFT/PolyPhen:** Variant effect predictions
- **Lange et al., 2005, Science:** Titin kinase signalosome
- **Bogomolovas et al., 2021, EMBO Rep:** TK ubiquitination and autophagy receptors

---

## License

Research project - AAAI-IAAA 2025

---

**Last Updated:** 2025-12-30
**Status:** ✅ Pipeline Complete (Steps 1-5) - R² > 0.99 on burden score prediction

