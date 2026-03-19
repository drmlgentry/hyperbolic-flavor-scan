import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ========= CONFIG =========
CSV_PATH = r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m006.csv"
OUT_DIR = Path(r"C:\dev\hyperbolic-flavor-scan\analysis")
BOOTSTRAPS = 1000
EPS = 1e-12
# ==========================

df = pd.read_csv(CSV_PATH)
OUT_DIR.mkdir(parents=True, exist_ok=True)

def running_min_curve(d):
    d = d.sort_values("length")
    d["phi_min"] = d["phi_fold"].cummin()
    g = d.groupby("length")["phi_min"].min().reset_index()
    return g["length"].to_numpy(), g["phi_min"].to_numpy()

def fit_exponential(L, phi):
    y = np.log(phi)
    b, a = np.polyfit(L, y, 1)
    c = -b
    A = np.exp(a)
    pred = A * np.exp(-c * L)
    return A, c, pred

def fit_powerlaw(L, phi):
    y = np.log(phi)
    x = np.log(L)
    b, a = np.polyfit(x, y, 1)
    alpha = -b
    A = np.exp(a)
    pred = A * L**(-alpha)
    return A, alpha, pred

def AIC(phi, pred, k):
    n = len(phi)
    rss = np.sum((phi - pred) ** 2)
    return n * np.log(rss / n) + 2 * k

summary = []

for k in sorted(df["class"].unique()):

    d = df[df["class"] == k].copy()
    L, phi = running_min_curve(d)

    mask = phi > EPS
    L, phi = L[mask], phi[mask]

    if len(L) < 5:
        continue

    # ---- Base fits ----
    Aexp, cexp, pred_exp = fit_exponential(L, phi)
    Apow, apow, pred_pow = fit_powerlaw(L, phi)

    aic_exp = AIC(phi, pred_exp, 2)
    aic_pow = AIC(phi, pred_pow, 2)

    # ---- Bootstrap ----
    boot_preds = []

    n = len(L)

    for _ in range(BOOTSTRAPS):
        idx = np.random.randint(0, n, n)
        Lb, phib = L[idx], phi[idx]

        order = np.argsort(Lb)
        Lb, phib = Lb[order], phib[order]

        try:
            A, c, pred = fit_exponential(Lb, phib)
            boot_preds.append(A * np.exp(-c * L))
        except:
            pass

    boot_preds = np.array(boot_preds)

    lo = np.percentile(boot_preds, 2.5, axis=0)
    hi = np.percentile(boot_preds, 97.5, axis=0)

    # ---- Plot ----
    plt.figure(figsize=(7, 5))

    plt.semilogy(L, phi, "o", label="data")
    plt.semilogy(L, pred_exp, "-", label="exp fit")
    plt.fill_between(L, lo, hi, alpha=0.25, label="95% CI")

    plt.xlabel("Word length L")
    plt.ylabel("Minimum twist (degrees)")
    plt.title(f"Class {k} decay")
    plt.grid(True, which="both", ls=":")
    plt.legend()

    out_path = OUT_DIR / f"decay_class_{k}.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    summary.append(
        (k, cexp, apow, aic_exp, aic_pow)
    )

# ---- Print model comparison ----
print("\nModel comparison (lower AIC = better)\n")
print("Class | c_exp | alpha_pow | AIC_exp | AIC_pow | Preferred")

for k, cexp, apow, aic_exp, aic_pow in summary:

    pref = "Exponential" if aic_exp < aic_pow else "Power-law"

    print(
        f"{k:5d} | {cexp:6.4f} | {apow:9.4f} | "
        f"{aic_exp:8.2f} | {aic_pow:8.2f} | {pref}"
    )

print(f"\nFigures saved in: {OUT_DIR}")