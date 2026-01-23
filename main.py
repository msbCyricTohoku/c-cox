import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.datasets import load_rossi
from scipy.stats import norm

rossi = load_rossi()
df = rossi[['week', 'arrest', 'fin', 'age']]

cph = CoxPHFitter() 
cph.fit(df, duration_col='week', event_col='arrest')

summary = cph.summary

print(f"\n{'Variable':<10} {'Coef':<10} {'SE':<10} {'p-val':<10} {'HR':<10} {'95% CI':<20}")
print("-" * 84)

for cov in ['fin', 'age']:
    beta = cph.params_[cov]
    se = np.sqrt(cph.variance_matrix_.loc[cov, cov])
    hr = np.exp(beta)
    
    # Wald Test
    z_stat = beta / se
    p_val = 2 * (1 - norm.cdf(abs(z_stat)))
    
    # 95% CI
    ci_low = np.exp(beta - 1.96 * se)
    ci_high = np.exp(beta + 1.96 * se)
    
    print(f"{cov:<10} {beta:<10.4f} {se:<10.4f} {p_val:<10.4f} {hr:<10.4f} ({ci_low:.4f}, {ci_high:.4f})")

print("-" * 84)

print("\n--- Variance-Covariance Matrix (Inverse Hessian) ---")
print(cph.variance_matrix_.loc[['fin', 'age'], ['fin', 'age']].to_string())
