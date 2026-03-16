"""
paper8_fig1.py
Figure 1: A-factor twist angle spectrum phi vs word length
for m003 and m006, with SM target lines overlaid.
Output: fig_twist_spectrum.pdf (vector, JHEP quality)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ── Load deduplicated census data ─────────────────────────────────────────────
m003 = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\twist_census_m003.csv",
                   skipinitialspace=True)
m006 = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\twist_census_m006.csv",
                   skipinitialspace=True)

m003.columns = [c.strip() for c in m003.columns]
m006.columns = [c.strip() for c in m006.columns]
m003["word_len"] = m003["word"].str.strip().str.len()
m006["word_len"] = m006["word"].str.strip().str.len()

# Deduplicate by phi (round to 3 dp), keep shortest word
def dedup(df):
    df = df.copy()
    df["phi_key"] = df["phi_deg"].round(3)
    df = df.sort_values(["phi_key", "word_len", "word"])
    return df.drop_duplicates(subset=["phi_key"], keep="first")

m003d = dedup(m003)
m006d = dedup(m006)

# Fold to [0, 90]: use min(|phi|, 180-|phi|) for display
def fold90(phi):
    p = abs(phi) % 180
    return min(p, 180 - p)

m003d = m003d.copy()
m006d = m006d.copy()
m003d["phi_fold"] = m003d["phi_deg"].apply(fold90)
m006d["phi_fold"] = m006d["phi_deg"].apply(fold90)

# ── SM target lines ───────────────────────────────────────────────────────────
SM_TARGETS = [
    (2.38,   r"$\theta_{23}^{\rm CKM}$",  "CKM",  "#1f77b4"),
    (13.04,  r"$\theta_{12}^{\rm CKM}$",  "CKM",  "#1f77b4"),
    (33.41,  r"$\theta_{12}^{\nu}$",       "PMNS", "#d62728"),
    (49.1,   r"$\theta_{23}^{\nu}$",       "PMNS", "#d62728"),
    (68.0,   r"$\delta_{\rm CKM}$",        "CKM",  "#1f77b4"),
    (8.54,   r"$\theta_{13}^{\nu}$",       "PMNS", "#d62728"),
]

# ── Figure setup ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
fig.subplots_adjust(wspace=0.08)

ALPHA_TARGET = 0.35
JITTER = 0.12   # horizontal jitter to separate overlapping points

for ax, df, title, color in [
    (axes[0], m003d, r"$M_{\rm PMNS} = \mathtt{m003}$", "#2ca02c"),
    (axes[1], m006d, r"$M_{\rm CKM}\ = \mathtt{m006}$",  "#ff7f0e"),
]:
    # SM target horizontal lines
    for tval, tlabel, sector, tcolor in SM_TARGETS:
        ax.axhline(tval, color=tcolor, alpha=ALPHA_TARGET,
                   linewidth=1.0, linestyle="--", zorder=1)
        ax.text(5.55, tval, tlabel, va="center", ha="left",
                fontsize=7, color=tcolor, alpha=0.85)

    # Scatter: jitter x slightly so stacked points separate
    rng = np.random.default_rng(42)
    x = df["word_len"] + rng.uniform(-JITTER, JITTER, len(df))
    ax.scatter(x, df["phi_fold"], color=color, s=22, zorder=3,
               edgecolors="white", linewidths=0.3, alpha=0.9)

    # Annotate the key hits
    KEY_HITS = {
        "m003": {
            "bbbb":  (r"$\mathtt{bbbb}$",    "mb/mc",        -3.5),
            "aabb":  (r"$\mathtt{aabb}$",    r"$\theta_{23}^\nu$",  2.5),
            "ABaB":  (r"$\mathtt{ABaB}$",    r"$\theta_{12}^\nu$",  2.5),
            "AAB":   (r"$\mathtt{AAB}$",     r"$\theta_{12}^{\rm CKM}$", -3.5),
        },
        "m006": {
            "aa":    (r"$\mathtt{aa}$",      r"$\delta_{\rm CKM}$",  2.5),
            "abbAB": (r"$\mathtt{abbAB}$",   r"$\theta_{12}^\nu$",  -3.5),
            "B":     (r"$\mathtt{b}$",       r"$\approx\pi/2$",      2.5),
            "aaabb": (r"$\mathtt{aaabb}$",   r"$\theta_{23}^{\rm CKM}$", 2.5),
        },
    }
    mname = "m003" if "m003" in title else "m006"
    for word_key, (wlabel, slabel, dy) in KEY_HITS.get(mname, {}).items():
        row = df[df["word"].str.strip().str.upper() == word_key.upper()]
        if row.empty:
            # try lowercase
            row = df[df["word"].str.strip() == word_key]
        if row.empty:
            continue
        xp = row["word_len"].values[0]
        yp = row["phi_fold"].values[0]
        ax.annotate(wlabel, xy=(xp, yp),
                    xytext=(xp + 0.25, yp + dy),
                    fontsize=6.5, color=color,
                    arrowprops=dict(arrowstyle="-", color=color,
                                   lw=0.6, alpha=0.7),
                    va="center")

    ax.set_title(title, fontsize=11, pad=6)
    ax.set_xlabel("Word length $|w|$", fontsize=10)
    ax.set_xlim(0.5, 5.7)
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_ylim(-1, 92)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, alpha=0.2, linewidth=0.5)

axes[0].set_ylabel(r"Folded twist angle $|\phi|$ (deg)", fontsize=10)

# Legend
ckm_line  = mlines.Line2D([], [], color="#1f77b4", linestyle="--",
                           alpha=0.7, label="CKM target")
pmns_line = mlines.Line2D([], [], color="#d62728", linestyle="--",
                           alpha=0.7, label="PMNS target")
fig.legend(handles=[ckm_line, pmns_line], loc="lower center",
           ncol=2, fontsize=9, framealpha=0.8,
           bbox_to_anchor=(0.5, -0.02))

fig.suptitle(r"A-factor twist angle spectrum: $\phi(\gamma)$ vs.\ word length",
             fontsize=12, y=1.01)

out = r"C:\dev\framework\papers\hyperbolic-flavor-twist\fig_twist_spectrum.pdf"
fig.savefig(out, format="pdf", bbox_inches="tight", dpi=300)
print(f"Saved: {out}")
