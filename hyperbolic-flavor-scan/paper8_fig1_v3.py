"""
paper8_fig1_v3.py
Figure 1 v3: broken axis, linear both panels, top-4 annotations.
Top panel: 25-92 deg. Bottom panel: 0-20 deg.
No log scale.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec

def load_dedup(path, name):
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    df["word"] = df["word"].str.strip()
    df["word_len"] = df["word"].str.len()
    df["phi_key"] = df["phi_deg"].round(3)
    df = df.sort_values(["phi_key", "word_len", "word"])
    df = df.drop_duplicates(subset=["phi_key"], keep="first")
    df["phi_fold"] = df["phi_deg"].apply(lambda p: min(abs(p)%180, 180-abs(p)%180))
    df["manifold"] = name
    return df

m003 = load_dedup(r"C:\dev\hyperbolic-flavor-scan\twist_census_m003.csv", "m003")
m006 = load_dedup(r"C:\dev\hyperbolic-flavor-scan\twist_census_m006.csv", "m006")

SM = [
    (2.38,  r"$\theta_{23}^{\rm CKM}$",  "#1f77b4"),
    (8.54,  r"$\theta_{13}^{\nu}$",       "#d62728"),
    (13.04, r"$\theta_{12}^{\rm CKM}$",  "#1f77b4"),
    (33.41, r"$\theta_{12}^{\nu}$",       "#d62728"),
    (49.1,  r"$\theta_{23}^{\nu}$",       "#d62728"),
    (68.0,  r"$\delta_{\rm CKM}$",        "#1f77b4"),
]

ANNOT = {
    "m003": [
        ("bbbb",  r"$\mathtt{bbbb}$",  "right",  5,   4),
        ("aabb",  r"$\mathtt{aabb}$",  "right",  5,   4),
        ("ABaB",  r"$\mathtt{ABaB}$",  "right",  5,  -5),
        ("AAB",   r"$\mathtt{AAB}$",   "right",  5,  -5),
    ],
    "m006": [
        ("aa",    r"$\mathtt{aa}$",    "left",  -5,   4),
        ("abbAB", r"$\mathtt{abbAB}$", "left",  -5,  -5),
        ("B",     r"$\mathtt{b}$",     "left",  -5,   4),
        ("aaabb", r"$\mathtt{aaabb}$", "left",  -5,  -5),
    ],
}

BREAK_LO = 21   # top of bottom panel
BREAK_HI = 27   # bottom of top panel

fig = plt.figure(figsize=(10, 6.5))
gs = GridSpec(2, 2, height_ratios=[2.6, 1.0], hspace=0.06, wspace=0.08,
              left=0.08, right=0.88, top=0.93, bottom=0.10)

axes_top = [fig.add_subplot(gs[0, i]) for i in range(2)]
axes_bot = [fig.add_subplot(gs[1, i]) for i in range(2)]

DATASETS = [
    (m003, "m003", r"$M_{\rm PMNS}=\mathtt{m003}$", "#2ca02c"),
    (m006, "m006", r"$M_{\rm CKM}\ =\mathtt{m006}$",  "#ff7f0e"),
]

for col, (df, mname, title, color) in enumerate(DATASETS):
    ax_t = axes_top[col]
    ax_b = axes_bot[col]

    rng = np.random.default_rng(7 + col)
    jitter = rng.uniform(-0.09, 0.09, len(df))
    x = df["word_len"] + jitter
    y = df["phi_fold"]

    # ── scatter both panels ────────────────────────────────────────────────
    ax_t.scatter(x, y, color=color, s=24, zorder=3,
                 edgecolors="white", linewidths=0.3, alpha=0.92)
    ax_b.scatter(x, y, color=color, s=24, zorder=3,
                 edgecolors="white", linewidths=0.3, alpha=0.92)

    # ── axis limits ────────────────────────────────────────────────────────
    ax_t.set_ylim(BREAK_HI, 92)
    ax_b.set_ylim(-0.5, BREAK_LO)

    for ax in [ax_t, ax_b]:
        ax.set_xlim(0.5, 5.55)
        ax.set_xticks([1, 2, 3, 4, 5])
        ax.grid(True, alpha=0.18, linewidth=0.5)

    ax_t.tick_params(bottom=False, labelbottom=False, labelsize=8)
    ax_b.tick_params(labelsize=8)
    ax_b.set_xlabel("Word length $|w|$", fontsize=10)
    ax_t.set_title(title, fontsize=11, pad=5)

    # ── SM target lines ────────────────────────────────────────────────────
    for tval, tlabel, tc in SM:
        for ax, lo, hi in [(ax_t, BREAK_HI, 92), (ax_b, -0.5, BREAK_LO)]:
            if lo <= tval <= hi:
                ax.axhline(tval, color=tc, alpha=0.38, lw=1.1, ls="--", zorder=1)
                ax.text(5.58, tval, tlabel, va="center", ha="left",
                        fontsize=7.5, color=tc, alpha=0.88)

    # ── broken axis marks ──────────────────────────────────────────────────
    ax_t.spines["bottom"].set_visible(False)
    ax_b.spines["top"].set_visible(False)
    d = 0.013
    kw = dict(color="k", clip_on=False, lw=0.9, transform=ax_t.transAxes)
    ax_t.plot((-d, +d), (-d*1.5, +d*1.5), **kw)
    ax_t.plot((1-d, 1+d), (-d*1.5, +d*1.5), **kw)
    kw.update(transform=ax_b.transAxes)
    ax_b.plot((-d, +d), (1-d*1.5, 1+d*1.5), **kw)
    ax_b.plot((1-d, 1+d), (1-d*1.5, 1+d*1.5), **kw)

    # ── annotations ────────────────────────────────────────────────────────
    for word_key, label, side, dx, dy in ANNOT.get(mname, []):
        row = df[df["word"].str.upper() == word_key.upper()]
        if row.empty:
            row = df[df["word"] == word_key]
        if row.empty:
            continue
        xp = float(row["word_len"].values[0])
        yp = float(row["phi_fold"].values[0])
        target_ax = ax_t if yp > BREAK_HI else ax_b
        target_ax.annotate(
            label,
            xy=(xp, yp),
            xytext=(xp + dx*0.12, yp + dy),
            fontsize=7, color=color,
            arrowprops=dict(arrowstyle="-", color=color, lw=0.55, alpha=0.7),
            va="center", ha="right" if side == "right" else "left",
        )

# ── shared y-axis labels ──────────────────────────────────────────────────────
fig.text(0.005, 0.70, r"$|\phi(\gamma)|$ (deg)",
         va="center", rotation="vertical", fontsize=10)
fig.text(0.005, 0.23, r"$|\phi(\gamma)|$ (deg)",
         va="center", rotation="vertical", fontsize=10)

# ── legend ────────────────────────────────────────────────────────────────────
ckm_l  = mlines.Line2D([], [], color="#1f77b4", ls="--", alpha=0.75,
                        label="CKM target (PDG 2024)")
pmns_l = mlines.Line2D([], [], color="#d62728", ls="--", alpha=0.75,
                        label="PMNS target (PDG 2024)")
fig.legend(handles=[ckm_l, pmns_l], loc="lower center", ncol=2,
           fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, 0.01))

fig.suptitle(
    r"A-factor twist spectrum: $|\phi(\gamma)|$ vs.\ word length",
    fontsize=12, y=0.98)

out = r"C:\dev\framework\papers\hyperbolic-flavor-twist\fig_twist_spectrum.pdf"
fig.savefig(out, format="pdf", bbox_inches="tight", dpi=300)
print(f"Saved: {out}")
