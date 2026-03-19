import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# ========= CONFIG =========
CSV_PATH = r"C:\dev\hyperbolic-flavor-scan\data\twist_random_L80_m006.csv"

OUT_PNG = r"C:\dev\hyperbolic-flavor-scan\analysis\twist_random_L80_hist.png"
EPS = 1e-12
# ==========================

df = pd.read_csv(CSV_PATH)

COL_CLASS = "h1_class"
COL_PHI = "phi_fold_deg"

print("\n===== RANDOM L=80 ANALYSIS =====\n")

# ---- Overall minimum ----
overall_min = df[COL_PHI].min()
print(f"Overall minimum twist: {overall_min:.6e} degrees")

# ---- Per-class minima ----
print("\nPer-class minima:")
mins = df.groupby(COL_CLASS)[COL_PHI].min()

for k, v in mins.items():
    print(f"Class {k}: {v:.6e} deg")

# ---- Histogram (log bins) ----
phi = df[COL_PHI].to_numpy()
phi = phi[phi > EPS]

bins = np.logspace(np.log10(phi.min()), np.log10(phi.max()), 60)

plt.figure(figsize=(8, 5))
plt.hist(phi, bins=bins)

plt.xscale("log")
plt.xlabel("Twist angle (degrees)")
plt.ylabel("Count")
plt.title("m006 — Random L=80 twist distribution")

Path(OUT_PNG).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")

plt.close("all")

print(f"\nSaved histogram: {OUT_PNG}")

sys.exit(0)