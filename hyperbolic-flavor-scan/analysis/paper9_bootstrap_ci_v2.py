import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

CSV_PATH = r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m006.csv"
OUT_DIR  = Path(r"C:\dev\hyperbolic-flavor-scan\analysis")
BOOTSTRAPS = 2000
EPS = 1e-12

df = pd.read_csv(CSV_PATH)
# Add translation length
df["trans_length"] = 2 * np.log(df["abs_lambda"].clip(lower=1.001))

def running_min_curve(d, xcol):
    d = d.sort_values(xcol)
    d = d.copy()
    d["phi_min"] = d["phi_fold_deg"].cummin()
    g = d.groupby(xcol)["phi_min"].min().reset_index()
    return g[xcol].to_numpy(), g["phi_min"].to_numpy()

def fit_exp(L, phi):
    y = np.log(phi + 1e-15)
    b, a = np.polyfit(L, y, 1)
    return np.exp(a), -b  # A, c

def AIC(phi, pred, k=2):
    n = len(phi)
    rss = np.sum((phi - pred)**2)
    return n * np.log(rss/n + 1e-30) + 2*k

print("="*70)
print("Word-length basis: c_k with 95% bootstrap CIs")
print("="*70)
print(f"{'Class':>6}  {'c_k':>8}  {'95% CI':>22}  {'ΔAIC':>8}  {'converged':>10}")

for basis, xcol in [("Word length", "length"), ("Translation length", "trans_length")]:
    print(f"\n--- {basis} basis ---")
    for k in range(5):
        d = df[df.h1_class == k].copy()
        L, phi = running_min_curve(d, xcol)
        mask = phi > EPS
        L, phi = L[mask], phi[mask]
        if len(L) < 4: continue

        A_hat, c_hat = fit_exp(L, phi)
        pred_exp = A_hat * np.exp(-c_hat * L)

        # Power-law for AIC comparison
        log_L = np.log(L + 1e-10)
        b, a = np.polyfit(log_L, np.log(phi+1e-15), 1)
        pred_pow = np.exp(a) * L**b
        delta_aic = AIC(phi, pred_exp) - AIC(phi, pred_pow)

        # Bootstrap c_k distribution
        rng = np.random.default_rng(42 + k)
        c_boots = []
        n = len(L)
        for _ in range(BOOTSTRAPS):
            idx = rng.integers(0, n, n)
            Lb, phib = L[idx], phi[idx]
            try:
                _, c_b = fit_exp(np.sort(Lb), phib[np.argsort(Lb)])
                if 0 < c_b < 20: c_boots.append(c_b)
            except: pass

        if len(c_boots) < 200:
            print(f"  Class {k}: bootstrap unstable ({len(c_boots)} converged)")
            continue

        lo, hi = np.percentile(c_boots, [2.5, 97.5])
        print(f"  Class {k}: c={c_hat:.3f}  [{lo:.3f}, {hi:.3f}]  "
              f"ΔAIC={delta_aic:+.1f}  ({len(c_boots)}/{BOOTSTRAPS})")

# Pairing verification
print("\n--- Pairing check ---")
print("c_1 and c_4 CIs should overlap (paired by k <-> 5-k involution)")
print("c_2 and c_3 CIs should overlap")
print("c_0 should be distinct from all others")
