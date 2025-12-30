# MechanobioAI: Methods, Results, and Discussion

**Supplementary Material for AAAI-IAAA 2025**

---

## Methods

### 1. Variant Data Extraction

Variant data for five proteins in the titin kinase signalosome (TTN, NBR1, p62/SQSTM1, MuRF1/TRIM63, MuRF2/TRIM55) were extracted from UniProt. Variants were filtered to those localized within annotated protein domains, yielding 2,139 domain-localized variants across 21 domain-protein combinations.

### 2. Impact Score Calculation

Each variant was assigned an impact score ∈ [0, 1] using a hierarchical approach:

| Source | Count | % | Formula |
|--------|-------|---|---------|
| Both SIFT & PolyPhen | 1,021 | 47.7% | `(1 - SIFT + PolyPhen) / 2` |
| SIFT only | 474 | 22.2% | `1 - SIFT` |
| PolyPhen only | 45 | 2.1% | `PolyPhen` |
| **Imputed** | **599** | **28.0%** | **Beta distribution by consequence type** |

For variants lacking both SIFT and PolyPhen scores, impact was sampled from Beta distributions parameterized by mutation consequence:
- **Stop-gained, frameshift, initiator codon**: Beta(8, 2), mean ≈ 0.80 (high impact)
- **Inframe deletion, delins**: Beta(6, 3), mean ≈ 0.67 (moderate-high)
- **Missense, insertion**: Beta(3, 3), mean ≈ 0.50 (variable)
- **Unknown**: Beta(2, 2), mean ≈ 0.50

### 3. Allele Frequency Imputation

Allele frequencies were obtained from gnomAD where available. For missing frequencies, values were sampled from a log-normal distribution fitted to observed rare-variant patterns:

```
log₁₀(frequency) ~ N(μ = -5.217, σ = 0.924)
```

This corresponds to a median frequency of ~6 × 10⁻⁶, consistent with rare variant expectations. Sampled values were clipped to the observed minimum (6.84 × 10⁻⁷) and maximum (1.0).

### 4. Domain-Pathway Weight Assignment

Weights linking each of 21 protein domains to hypertrophy and autophagy pathways were assigned via LLM-assisted literature review. Two LLMs (Claude Opus 4 and GPT-4o) independently assigned weights ∈ [0, 1] based on domain function and primary literature references (Lange et al. 2005, Bogomolovas et al. 2021, Bodine et al. 2001, etc.). Final weights were computed as the average of both models. Each weight assignment includes two literature references with relevance descriptions.

Example weights:
- **TTN Kinase → Hypertrophy**: 0.80 (central mechanosensing hub)
- **TTN Kinase → Autophagy**: 0.825 (recruits MuRF1/NBR1/p62 complex)
- **p62 LIR → Autophagy**: 0.925 (essential autophagosome targeting)
- **MuRF1 RING → Hypertrophy**: 0.80 (master atrophy regulator)

### 5. Synthetic Cohort Generation

For n=50,000 synthetic individuals, genotypes were sampled as:
```
g_i ~ Binomial(2, allele_freq_i)  for each variant i
```

Pathway burden scores were computed as:
```
hypertrophy_burden = Σᵢ (impact_i × hyp_weight_i × g_i)
autophagy_burden = Σᵢ (impact_i × auto_weight_i × g_i)
```

### 6. Phenotype Generation

Ten observable phenotypes were generated as linear combinations of pathway burdens:
```
phenotype = baseline + β_hyp × hyp_burden + β_auto × auto_burden + ε
```
where ε ~ N(0, σ²) with phenotype-specific noise.

| Phenotype | β_hyp | β_auto | σ | Rationale |
|-----------|-------|--------|---|-----------|
| muscle_growth_rate | -1.5 | +0.5 | 0.12 | Hyp disruption impairs growth; auto disruption reduces catabolism |
| natural_muscle_mass | -1.5 | +0.6 | 0.10 | Same logic; autophagy disruption preserves mass |
| strength_1rm | -1.2 | +0.3 | 0.15 | Strength correlates with muscle mass |
| power_output | -1.3 | +0.2 | 0.12 | Fast-twitch requires hypertrophic signaling |
| vo2max | -0.2 | -1.5 | 0.15 | Autophagy critical for mitochondrial quality |
| time_to_exhaustion | -0.1 | -1.5 | 0.12 | Mitophagy maintains aerobic capacity |
| recovery_speed | -0.2 | -1.5 | 0.12 | Autophagy clears exercise-induced damage |
| doms_severity | +0.2 | +1.5 | 0.15 | Impaired autophagy = more soreness |
| injury_susceptibility | +0.4 | +0.8 | 0.15 | Both pathways affect tissue quality |
| training_responsiveness | -0.6 | -0.6 | 0.12 | Both equally required for adaptation |

### 7. Inverse Prediction Model

A multi-layer perceptron was trained to predict the two latent burden scores from the 10 observable phenotypes:

- **Architecture**: 10 → 64 → 32 → 2 (ReLU activation)
- **Optimizer**: Adam with early stopping (patience=20, validation=10%)
- **Split**: 85% train (42,500), 15% test (7,500)
- **Preprocessing**: StandardScaler on input features
- **Seed**: 42 for full reproducibility

---

## Results

The trained MLP achieved near-perfect recovery of latent pathway burdens from observable phenotypes, with test R² = 0.990 for hypertrophy burden and R² = 0.994 for autophagy burden (n=7,500 test samples). Train-test performance was virtually identical, indicating no overfitting. Mean absolute prediction errors were 0.019 and 0.020 for hypertrophy and autophagy burdens respectively, on burden scales with standard deviations of 0.31 and 0.39.

Partial correlation analysis confirmed pathway-specific phenotype associations: muscle phenotypes (growth rate, mass, strength, power) showed strong unique correlations with hypertrophy burden (|r| = 0.75–0.91) but weaker or opposite correlations with autophagy. Conversely, endurance and recovery phenotypes (VO2max, time-to-exhaustion, recovery speed, DOMS) were dominated by autophagy burden (|r| = 0.87–0.91). Training responsiveness showed balanced contributions from both pathways. These patterns align with known biology: hypertrophy pathway disruption impairs anabolic muscle-building, while autophagy disruption impairs mitochondrial quality control critical for endurance and recovery.

---

## Discussion

This feasibility study demonstrates that mechanobiology-informed phenotype patterns can, in principle, carry recoverable signal about underlying pathway-level genetic burden. However, several critical limitations prevent direct clinical translation.

**Biophysical modeling gap:** Our impact scores derive from sequence-based predictors (SIFT/PolyPhen), not from actual biophysical simulations of how mutations alter protein dynamics or binding affinities. True mechanobiological modeling would require structure prediction for mutated proteins (challenging for giant proteins like titin) followed by molecular dynamics simulations to quantify functional effects on mechanosensing.

**Knowledge graph incompleteness:** Domain-pathway and genomic trait-phenotype weights were assigned via LLM-assisted literature review rather than systematic ontology extraction. A production system would require comprehensive mining of PubMed to build exhaustive protein-pathway-phenotype relationships with quantified confidence.

**Absence of real validation data:** No datasets currently pair individual genotypes with controlled exercise-adaptation outcomes at sufficient scale. The synthetic cohort validates computational consistency but cannot confirm biological validity. Prospective studies pairing genomic profiles with longitudinal training responses are essential for model validation.

Despite these limitations, the framework establishes a principled architecture for N-of-1 exercise personalization: observable phenotypes → inferred pathway burdens → targeted loading prescriptions designed to activate specific mechanotransduction cascades. Addressing the above gaps would transform this proof-of-concept into a clinically actionable system.

---

## References

1. Ibata & Terentjev (2021). Biophysical Journal. Titin mechanosensing and muscle growth.
2. Lange et al. (2005). Science. Titin kinase signalosome.
3. Bogomolovas et al. (2021). EMBO Reports. TK ubiquitination and autophagy receptors.
4. Bodine et al. (2001). Science. MuRF1 and muscle atrophy.
5. Sandri et al. (2013). Biogerontology. Autophagy in muscle homeostasis.

