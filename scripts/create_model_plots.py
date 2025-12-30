#!/usr/bin/env python3
"""
Create model training and evaluation plots.
Run this AFTER train_burden_predictor.py
"""

import numpy as np
import matplotlib.pyplot as plt
import json

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150

print("Creating model evaluation plots...")

# Load saved data
loss_curve = np.load('../synthetic_cohort/loss_curve.npy')
y_test = np.load('../synthetic_cohort/y_test.npy')
y_test_pred = np.load('../synthetic_cohort/y_test_pred.npy')
y_train = np.load('../synthetic_cohort/y_train.npy')
y_train_pred = np.load('../synthetic_cohort/y_train_pred.npy')

with open('../synthetic_cohort/model_results.json', 'r') as f:
    results = json.load(f)

target_cols = ['hypertrophy_burden', 'autophagy_burden']

# =============================================================================
# PLOT 1: Training Loss Curve
# =============================================================================
fig1, ax1 = plt.subplots(figsize=(10, 5))

ax1.plot(loss_curve, color='#1976D2', linewidth=2, label='Training Loss')
ax1.set_xlabel('Epoch')
ax1.set_ylabel('Loss (MSE)')
ax1.set_title('Training Loss Curve\n(MLP: 64-32 hidden units)', fontsize=12, fontweight='bold')
ax1.set_yscale('log')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Add annotation
n_iter = results['model_info']['n_iter']
final_loss = loss_curve[-1]
ax1.annotate(f'Final: {final_loss:.4f}\n({n_iter} epochs)', 
             xy=(n_iter-1, final_loss), xytext=(n_iter*0.7, final_loss*2),
             arrowprops=dict(arrowstyle='->', color='gray'),
             fontsize=10)

fig1.tight_layout()
fig1.savefig('../plots/model_training_loss.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved model_training_loss.png")
plt.close(fig1)

# =============================================================================
# PLOT 2: Predicted vs True (Test Set)
# =============================================================================
fig2, axes = plt.subplots(1, 2, figsize=(12, 5))

for i, (ax, target) in enumerate(zip(axes, target_cols)):
    true = y_test[:, i]
    pred = y_test_pred[:, i]
    r2 = results[target]['r2_test']
    
    ax.scatter(true, pred, alpha=0.3, s=5, c='#1976D2')
    
    # Perfect prediction line
    lims = [min(true.min(), pred.min()), max(true.max(), pred.max())]
    ax.plot(lims, lims, 'r--', linewidth=2, label='Perfect prediction')
    
    ax.set_xlabel(f'True {target.replace("_", " ").title()}')
    ax.set_ylabel(f'Predicted {target.replace("_", " ").title()}')
    ax.set_title(f'{target.replace("_", " ").title()}\nTest R² = {r2:.4f}')
    ax.legend(loc='lower right')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

fig2.suptitle('Predicted vs True Burden Scores (Test Set, n=7,500)', fontsize=12, fontweight='bold')
fig2.tight_layout()
fig2.savefig('../plots/model_pred_vs_true.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved model_pred_vs_true.png")
plt.close(fig2)

# =============================================================================
# PLOT 3: Error Distribution
# =============================================================================
fig3, axes = plt.subplots(1, 2, figsize=(12, 5))

for i, (ax, target) in enumerate(zip(axes, target_cols)):
    errors = y_test_pred[:, i] - y_test[:, i]
    mae = results[target]['mae_test']
    
    ax.hist(errors, bins=50, color='#388E3C', alpha=0.7, edgecolor='white')
    ax.axvline(0, color='red', linestyle='--', linewidth=2)
    ax.axvline(errors.mean(), color='blue', linestyle='-', linewidth=2, label=f'Mean: {errors.mean():.4f}')
    
    ax.set_xlabel('Prediction Error')
    ax.set_ylabel('Count')
    ax.set_title(f'{target.replace("_", " ").title()}\nMAE = {mae:.4f}')
    ax.legend()

fig3.suptitle('Prediction Error Distribution (Test Set)', fontsize=12, fontweight='bold')
fig3.tight_layout()
fig3.savefig('../plots/model_error_distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved model_error_distribution.png")
plt.close(fig3)

# =============================================================================
# PLOT 4: Train vs Test Performance Summary
# =============================================================================
fig4, ax4 = plt.subplots(figsize=(10, 6))

metrics = ['R²', 'MSE', 'MAE']
x = np.arange(len(target_cols))
width = 0.35

# R² comparison
r2_train = [results[t]['r2_train'] for t in target_cols]
r2_test = [results[t]['r2_test'] for t in target_cols]

bars1 = ax4.bar(x - width/2, r2_train, width, label='Train', color='#1976D2', alpha=0.8)
bars2 = ax4.bar(x + width/2, r2_test, width, label='Test', color='#D32F2F', alpha=0.8)

ax4.set_ylabel('R² Score')
ax4.set_title('Model Performance: Train vs Test\n(Higher is better)', fontsize=12, fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels([t.replace('_', '\n').title() for t in target_cols])
ax4.legend()
ax4.set_ylim(0, 1)

# Add value labels
for bar, val in zip(bars1, r2_train):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{val:.3f}', 
             ha='center', va='bottom', fontsize=10)
for bar, val in zip(bars2, r2_test):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{val:.3f}',
             ha='center', va='bottom', fontsize=10)

fig4.tight_layout()
fig4.savefig('../plots/model_performance_summary.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  Saved model_performance_summary.png")
plt.close(fig4)

print("\nModel evaluation plots complete!")
print(f"\nSummary:")
for target in target_cols:
    print(f"  {target}: Train R²={results[target]['r2_train']:.3f}, Test R²={results[target]['r2_test']:.3f}")

