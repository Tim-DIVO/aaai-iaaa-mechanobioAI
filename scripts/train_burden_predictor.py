#!/usr/bin/env python3
"""
Train a simple MLP to predict burden scores from observable phenotypes.
- Input: 10 phenotypes
- Output: 2 burden scores (hypertrophy, autophagy)
- Model: 2 hidden layers, no hyperparameter tuning
- Split: 85/15 train/test
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import json

# Config
SEED = 42
TEST_SIZE = 0.15
HIDDEN_LAYERS = (64, 32)  # 2 hidden layers
MAX_ITER = 500

np.random.seed(SEED)

# Load data
print("Loading synthetic cohort...")
cohort = pd.read_csv('../synthetic_cohort/synthetic_cohort_summary.csv')
phenotype_defs = pd.read_csv('../trait_mapping/phenotype_definitions.csv')

# Define features and targets
phenotype_cols = phenotype_defs['phenotype'].tolist()
target_cols = ['hypertrophy_burden', 'autophagy_burden']

X = cohort[phenotype_cols].values
y = cohort[target_cols].values

print(f"Features: {len(phenotype_cols)} phenotypes")
print(f"Targets: {len(target_cols)} burden scores")
print(f"Samples: {len(X)}")

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=TEST_SIZE, random_state=SEED
)
print(f"\nTrain: {len(X_train)}, Test: {len(X_test)}")

# Standardize features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train MLP
print(f"\nTraining MLP with hidden layers: {HIDDEN_LAYERS}...")
model = MLPRegressor(
    hidden_layer_sizes=HIDDEN_LAYERS,
    activation='relu',
    solver='adam',
    max_iter=MAX_ITER,
    random_state=SEED,
    early_stopping=True,
    validation_fraction=0.1,
    n_iter_no_change=20,
    verbose=True
)

model.fit(X_train_scaled, y_train)

# Predictions
y_train_pred = model.predict(X_train_scaled)
y_test_pred = model.predict(X_test_scaled)

# Evaluate
print("\n" + "="*60)
print("MODEL EVALUATION")
print("="*60)

results = {}
for i, target in enumerate(target_cols):
    r2_train = r2_score(y_train[:, i], y_train_pred[:, i])
    r2_test = r2_score(y_test[:, i], y_test_pred[:, i])
    mse_train = mean_squared_error(y_train[:, i], y_train_pred[:, i])
    mse_test = mean_squared_error(y_test[:, i], y_test_pred[:, i])
    mae_train = mean_absolute_error(y_train[:, i], y_train_pred[:, i])
    mae_test = mean_absolute_error(y_test[:, i], y_test_pred[:, i])
    
    results[target] = {
        'r2_train': r2_train, 'r2_test': r2_test,
        'mse_train': mse_train, 'mse_test': mse_test,
        'mae_train': mae_train, 'mae_test': mae_test
    }
    
    print(f"\n{target}:")
    print(f"  R² (train): {r2_train:.4f}")
    print(f"  R² (test):  {r2_test:.4f}")
    print(f"  MSE (test): {mse_test:.6f}")
    print(f"  MAE (test): {mae_test:.4f}")

# Save results
results['model_info'] = {
    'hidden_layers': HIDDEN_LAYERS,
    'n_iter': model.n_iter_,
    'train_size': len(X_train),
    'test_size': len(X_test),
    'features': phenotype_cols,
    'targets': target_cols
}

with open('../synthetic_cohort/model_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print("\nSaved results to ../synthetic_cohort/model_results.json")

# Save loss curve
loss_curve = model.loss_curve_
val_scores = model.validation_scores_ if hasattr(model, 'validation_scores_') else None

np.save('../synthetic_cohort/loss_curve.npy', loss_curve)
if val_scores:
    np.save('../synthetic_cohort/val_scores.npy', val_scores)

# Save predictions for plotting
np.save('../synthetic_cohort/y_test.npy', y_test)
np.save('../synthetic_cohort/y_test_pred.npy', y_test_pred)
np.save('../synthetic_cohort/y_train.npy', y_train)
np.save('../synthetic_cohort/y_train_pred.npy', y_train_pred)

print("\nTraining complete!")
print(f"Epochs: {model.n_iter_}")
print(f"Final loss: {loss_curve[-1]:.6f}")

