import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

df = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m006.csv")
# columns: word, length, h1_class, phi_fold_deg, abs_lambda
df["trans_length"] = 2 * np.log(df["abs_lambda"].clip(lower=1.001))

COLORS = {0:"#1f77b4", 1:"#ff7f0e", 2:"#2ca02c", 3:"#d62728", 4:"#9467bd"}
LABELS = {0:"Class 0", 1:"Class 1", 2:"Class 2", 3:"Class 3", 4:"Class 4"}

fig, ax = plt.subplots(figsize=(8, 5.5))

for k in range(5):
    sub = df[df.h1_class == k].sort_values("trans_length")
    # cumulative floor vs max translation length seen so far
    floors, mids = [], []
    current_floor = np.inf
    for _, row in sub.iterrows():
        if row.phi_fold_deg < current_floor:
            current_floor = row.phi_fold_deg
        floors.append(current_floor)
        mids.append(row.trans_length)
    
    floors = np.array(floors)
    mids   = np.array(mids)
    mask   = floors > 1e-4
    
    # Downsample for plot clarity
    step = max(1, len(mids[mask])//80)
    ax.scatter(mids[mask][::step], floors[mask][::step],
               color=COLORS[k], s=20, alpha=0.75, label=LABELS[k], zorder=3)
    
    # Exponential fit on downsampled
    try:
        def exp_f(x, A, c): return A * np.exp(-c * x)
        popt, _ = curve_fit(exp_f, mids[mask][::step], floors[mask][::step],
                            p0=[50, 0.3], maxfev=5000)
        x_fit = np.linspace(mids[mask].min(), mids[mask].max(), 200)
        ax.plot(x_fit, exp_f(x_fit, *popt), color=COLORS[k],
                ls="--", lw=1.2, alpha=0.8)
        print(f"Class {k}: c={popt[1]:.3f} (translation length basis)")
    except Exception as e:
        print(f"Class {k}: fit failed: {e}")

ax.set_yscale("log")
ax.set_xlabel(r"Cumulative max translation length $\max_{\ell\leq L}\ell(\gamma)=2\log|\lambda|$",
              fontsize=10)
ax.set_ylabel(r"$\phi_{\rm floor}^{(k)}$ (deg)", fontsize=11)
ax.set_title(r"m006: spectral floor vs. translation length (generator-independent)",
             fontsize=11)
ax.grid(True, alpha=0.15, lw=0.5, which="both")
ax.legend(fontsize=9, loc="upper right")
fig.tight_layout()

out = r"C:\dev\hyperbolic-flavor-scan\figures\fig_floor_decay_translength.pdf"
fig.savefig(out, format="pdf", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved: {out}")
