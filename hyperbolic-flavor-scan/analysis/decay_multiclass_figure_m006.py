import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# ========= CONFIG =========
CSV_PATH = r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m006.csv"

OUT_PNG = r"C:\dev\hyperbolic-flavor-scan\analysis\decay_multiclass_m006_clean.png"
OUT_PDF = r"C:\dev\hyperbolic-flavor-scan\analysis\decay_multiclass_m006_clean.pdf"

SHOW_PLOT = False   # ← set True only if you want a quick preview
EPS = 1e-12
JITTER = 0.06
# ==========================

df = pd.read_csv(CSV_PATH)

COL_LEN = "length"
COL_CLASS = "h1_class"
COL_PHI = "phi_fold_deg"

colors = {
    0: "tab:blue",
    1: "#FFD700",   # gold
    2: "tab:green",
    3: "tab:purple",
    4: "#DC143C"    # scarlet
}

offsets = {1: -JITTER, 4: +JITTER}

fig, ax = plt.subplots(figsize=(8.5, 6.2))

for k in sorted(df[COL_CLASS].unique()):

    d = df[df[COL_CLASS] == k].copy()
    d = d.sort_values(COL_LEN)

    d["phi_min"] = d[COL_PHI].cummin()
    g = d.groupby(COL_LEN)["phi_min"].min().reset_index()

    L = g[COL_LEN].to_numpy()
    phi = g["phi_min"].to_numpy()

    mask = phi > EPS
    L, phi = L[mask], phi[mask]

    if len(L) < 3:
        continue

    L_plot = L + offsets.get(k, 0)

    # Exponential fit
    y = np.log(phi)
    b, a = np.polyfit(L, y, 1)
    c = -b
    A = np.exp(a)

    color = colors.get(k, "black")

    ax.semilogy(
        L_plot, phi, "o",
        color=color,
        markersize=7,
        markeredgecolor="black",
        markeredgewidth=0.6,
        label=f"class {k}"
    )

    Lfit = np.linspace(min(L), max(L), 300)
    ax.semilogy(
        Lfit,
        A * np.exp(-c * Lfit),
        "--",
        color=color,
        linewidth=1.6
    )

ax.set_xlabel("Word length $L$")
ax.set_ylabel("Minimum twist (degrees)")
ax.set_title("m006 — Exponential decay of twist floors by homology class")

ax.grid(True, which="both", linestyle=":", alpha=0.7)
ax.legend(title="Homology class", framealpha=0.95)

# ---- Save outputs ----
Path(OUT_PNG).parent.mkdir(parents=True, exist_ok=True)

fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_PDF, bbox_inches="tight")

# ---- Optional preview (non-blocking) ----
if SHOW_PLOT:
    plt.show(block=False)
    plt.pause(2)   # show for 2 seconds

# ---- Clean exit ----
plt.close('all')
sys.exit(0)