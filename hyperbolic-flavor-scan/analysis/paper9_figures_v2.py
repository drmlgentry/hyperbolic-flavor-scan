"""
paper9_figures_v2.py - Updated figures with length-10 filtered data
Fig 1: phi_fold vs word length, colored by homology class (m006, lengths 1-10)
Fig 2: Bar chart of spectral floors - length-7 (m003) and length-10 filtered (m006)
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches

def homology_class(word, n=5):
    return sum(1 if c.islower() else -1 for c in word) % n

COSET_COLORS = {0:"#1f77b4", 1:"#ff7f0e", 2:"#2ca02c", 3:"#d62728", 4:"#9467bd"}

# ── Load length-10 data (m006) ────────────────────────────────────
df10 = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len10_m006.csv")
df10["h1_class"] = df10["word"].apply(homology_class)

# Filter trivials
df_geo = df10[df10.mod_lambda > 1.01].copy()

# Deduplicate by phi_fold (keep shortest word per unique value)
df_geo["phi_key"] = df_geo["phi_fold"].round(3)
df_geo = df_geo.sort_values(["phi_key","length","word"])
df_geo_dedup = df_geo.drop_duplicates("phi_key", keep="first")

# ── Figure 1: scatter phi_fold vs length, colored by coset ────────
fig, ax = plt.subplots(figsize=(8, 5))
rng = np.random.default_rng(42)

for coset in range(5):
    sub = df_geo_dedup[df_geo_dedup.h1_class == coset]
    jitter = rng.uniform(-0.10, 0.10, len(sub))
    ax.scatter(sub["length"] + jitter, sub["phi_fold"],
               color=COSET_COLORS[coset], s=18, alpha=0.80,
               edgecolors="none", label=f"Class {coset}", zorder=3)

# Annotate key points
annots = [
    ("AAABAB",  4, "right", r"$\mathtt{AAABAB}$ (cls 4)"),
    ("aaabab",  1, "left",  r"$\mathtt{aaabab}$ (cls 1)"),
]
for word, cls, ha, label in annots:
    row = df_geo_dedup[df_geo_dedup.word == word]
    if len(row):
        xp, yp = float(row.length), float(row.phi_fold)
        ax.annotate(label, xy=(xp, yp),
                    xytext=(xp + (1.0 if ha=="left" else -1.0), yp + 3.5),
                    fontsize=7.5, color=COSET_COLORS[cls], ha=ha,
                    arrowprops=dict(arrowstyle="-", color=COSET_COLORS[cls], lw=0.6))

# Gap shading [0.005, 0.253]
ax.axhspan(0.005, 0.253, alpha=0.07, color="gray")
ax.text(10.15, 0.13, "Gap\n$[0.005°,0.253°]$", va="center", ha="left",
        fontsize=7, color="gray")

# SM target
ax.axhline(0.201, color="gray", ls="--", lw=0.9, alpha=0.7)
ax.text(10.15, 0.201, r"$\theta_{13}^{\rm CKM}$", va="center",
        ha="left", fontsize=8, color="gray")

ax.set_xlabel("Word length $|w|$", fontsize=11)
ax.set_ylabel(r"$|\phi(\gamma)|$ (deg)", fontsize=11)
ax.set_title(r"\texttt{m006}: twist spectrum by homology class ($|w|\leq 10$, genuine loxodromics)",
             fontsize=10)
ax.set_xlim(0.5, 10.7)
ax.set_xticks(range(1,11))
ax.set_ylim(-1, 92)
ax.grid(True, alpha=0.12, lw=0.5)
ax.legend(fontsize=8, loc="upper left", framealpha=0.85, ncol=3)

# Add note about trivials
ax.text(0.98, 0.02,
        "18 trivial elements ($|\\lambda|\\leq 1.01$) excluded",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=7, color="gray", style="italic")

fig.tight_layout()
fig.savefig(r"C:\dev\framework\papers\hyperbolic-flavor-torsion\fig_coset_spectrum.pdf",
            format="pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved fig_coset_spectrum.pdf")

# ── Figure 2: bar chart - m003 len7, m006 len10 filtered ──────────
floors_m003 = {0:2.643, 1:2.763, 2:4.857, 3:5.376, 4:3.269}
floors_m006 = {0:1.023, 1:0.005, 2:0.253, 3:0.253, 4:0.005}  # len-10 filtered

fig, axes = plt.subplots(1, 2, figsize=(9, 4))

data = [
    (axes[0], floors_m003, r"\texttt{m003} — Meyerhoff ($|w|\leq 7$)", "Max/min: $2\\times$"),
    (axes[1], floors_m006, r"\texttt{m006} ($|w|\leq 10$, genuine loxodromics)",
     "Gap $[0.005°, 0.253°]$\nMax/min: $193\\times$"),
]

for ax, floors, title, note in data:
    cosets = list(floors.keys())
    vals   = [floors[k] for k in cosets]
    bars   = ax.bar(cosets, vals,
                    color=[COSET_COLORS[k] for k in cosets],
                    edgecolor="white", linewidth=0.5, zorder=3)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + max(vals)*0.02,
                f"{val:.3f}°", ha="center", va="bottom", fontsize=7.5)
    ax.set_xlabel("Homology class $k\\in\\mathbb{Z}/5$", fontsize=10)
    ax.set_ylabel(r"$\phi_{\rm floor}(k)$ (deg)", fontsize=10)
    ax.set_title(title, fontsize=9.5)
    ax.set_xticks(range(5))
    ax.grid(True, alpha=0.15, lw=0.5, axis="y")
    ax.set_ylim(0, max(vals)*1.32)
    ax.text(0.97, 0.97, note, transform=ax.transAxes,
            ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                      edgecolor="gray", alpha=0.85))

# Shade the gap region on m006 panel
axes[1].axhspan(0.005, 0.253, alpha=0.10, color="gray", zorder=0)
axes[1].axhline(0.201, color="gray", ls="--", lw=0.9, alpha=0.7)
axes[1].text(4.6, 0.215, r"$\theta_{13}^{\rm CKM}$",
             fontsize=7.5, color="gray", ha="right")

fig.suptitle("Spectral floors per homology class", fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(r"C:\dev\framework\papers\hyperbolic-flavor-torsion\fig_coset_floors.pdf",
            format="pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved fig_coset_floors.pdf")
