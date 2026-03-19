"""
paper9_figures.py - Generate figures for Paper 9
Fig 1: phi_fold vs word length, colored by homology class (m006)
Fig 2: Bar chart of spectral floors per coset (both manifolds)
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def homology_class(word, n=5):
    return sum(1 if c.islower() else -1 for c in word) % n

COSET_COLORS = {0:"#1f77b4", 1:"#ff7f0e", 2:"#2ca02c", 3:"#d62728", 4:"#9467bd"}
COSET_LABELS = {k: f"Class {k}" for k in range(5)}

# ── Load m006 length-7 data ───────────────────────────────────────
df6 = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len7_m006.csv")
df6["phi_fold"] = df6["phi_deg"].apply(lambda p: min(abs(p)%180, 180-abs(p)%180))
df6["phi_key"]  = df6["phi_fold"].round(3)
df6 = df6.sort_values(["phi_key","length","word"]).drop_duplicates("phi_key", keep="first")
df6["h1_class"] = df6["word"].apply(homology_class)

# ── Figure 1: scatter phi_fold vs length, colored by coset ────────
fig, ax = plt.subplots(figsize=(7, 4.5))
rng = np.random.default_rng(42)

for coset in range(5):
    sub = df6[df6.h1_class == coset]
    jitter = rng.uniform(-0.08, 0.08, len(sub))
    ax.scatter(sub["length"] + jitter, sub["phi_fold"],
               color=COSET_COLORS[coset], s=28, alpha=0.85,
               edgecolors="white", linewidths=0.3,
               label=COSET_LABELS[coset], zorder=3)

# Annotate the isolated coset-4 point
low = df6[df6.phi_fold < 0.1]
for _, row in low.iterrows():
    ax.annotate(r"$\mathtt{AAABAB}$" + f"\n({row.phi_fold:.4f}°)",
                xy=(row.length, row.phi_fold),
                xytext=(row.length - 1.2, 4.0),
                fontsize=7.5, color=COSET_COLORS[4],
                arrowprops=dict(arrowstyle="-", color=COSET_COLORS[4], lw=0.7))

# SM target line
ax.axhline(0.201, color="gray", ls="--", lw=0.9, alpha=0.6)
ax.text(7.08, 0.201, r"$\theta_{13}^{\rm CKM}$", va="center",
        ha="left", fontsize=8, color="gray")

# Gap shading
ax.axhspan(0.005, 1.611, alpha=0.07, color="gray", label="Spectral gap")

ax.set_xlabel("Word length $|w|$", fontsize=11)
ax.set_ylabel(r"$|\phi(\gamma)|$ (deg)", fontsize=11)
ax.set_title(r"\texttt{m006}: twist spectrum by homology class", fontsize=11)
ax.set_xlim(0.5, 7.7)
ax.set_xticks(range(1,8))
ax.set_ylim(-1, 92)
ax.grid(True, alpha=0.15, lw=0.5)
ax.legend(fontsize=8, loc="upper left", framealpha=0.85, ncol=2)
fig.tight_layout()
fig.savefig(r"C:\dev\framework\papers\hyperbolic-flavor-torsion\fig_coset_spectrum.pdf",
            format="pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved fig_coset_spectrum.pdf")

# ── Figure 2: bar chart of spectral floors ────────────────────────
floors_m003 = {0:2.643, 1:2.763, 2:4.857, 3:5.376, 4:3.269}
floors_m006 = {0:2.132, 1:3.374, 2:2.399, 3:1.611, 4:0.005}

fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))

for ax, floors, name, color in [
    (axes[0], floors_m003, r"\texttt{m003} (Meyerhoff)", "#2ca02c"),
    (axes[1], floors_m006, r"\texttt{m006}", "#ff7f0e")]:

    cosets = list(floors.keys())
    vals   = [floors[k] for k in cosets]
    bars   = ax.bar(cosets, vals,
                    color=[COSET_COLORS[k] for k in cosets],
                    edgecolor="white", linewidth=0.5, zorder=3)

    # Label bars
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f"{val:.3f}°", ha="center", va="bottom", fontsize=7.5)

    ax.set_xlabel("Homology class $k \\in \\mathbb{Z}/5$", fontsize=10)
    ax.set_ylabel(r"$\phi_{\rm floor}(k)$ (deg)", fontsize=10)
    ax.set_title(name, fontsize=10)
    ax.set_xticks(range(5))
    ax.grid(True, alpha=0.15, lw=0.5, axis="y")
    ax.set_ylim(0, max(vals)*1.25)

# Add ratio annotation to m006 panel
axes[1].text(0.97, 0.97, "Max/min ratio: $637\\times$",
             transform=axes[1].transAxes, ha="right", va="top",
             fontsize=8.5, color="darkred",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow",
                       edgecolor="darkred", alpha=0.8))

fig.suptitle("Spectral floors per homology class", fontsize=11, y=1.01)
fig.tight_layout()
fig.savefig(r"C:\dev\framework\papers\hyperbolic-flavor-torsion\fig_coset_floors.pdf",
            format="pdf", dpi=300, bbox_inches="tight")
plt.close()
print("Saved fig_coset_floors.pdf")
