"""
paper8_fig1_v2.py
Figure 1 v2: log y-scale, top-4 annotations only per panel,
broken axis visual cue between 15 and 28 deg.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import ConnectionPatch

# ── Load and deduplicate ──────────────────────────────────────────────────────
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

# ── SM targets ────────────────────────────────────────────────────────────────
SM = [
    (2.38,  r"$\theta_{23}^{\rm CKM}$",  "#1f77b4"),
    (8.54,  r"$\theta_{13}^{\nu}$",       "#d62728"),
    (13.04, r"$\theta_{12}^{\rm CKM}$",  "#1f77b4"),
    (33.41, r"$\theta_{12}^{\nu}$",       "#d62728"),
    (49.1,  r"$\theta_{23}^{\nu}$",       "#d62728"),
    (68.0,  r"$\delta_{\rm CKM}$",        "#1f77b4"),
]

# Top-4 annotation targets per manifold (word, label, dy_points)
ANNOT = {
    "m003": [
        ("bbbb",  r"$\mathtt{bbbb}$: $m_b/m_c$",            10),
        ("aabb",  r"$\mathtt{aabb}$: $\theta_{23}^\nu$",     10),
        ("ABaB",  r"$\mathtt{ABaB}$: $\theta_{12}^\nu$",    -12),
        ("AAB",   r"$\mathtt{AAB}$: $\theta_{12}^{\rm CKM}$", -12),
    ],
    "m006": [
        ("aa",    r"$\mathtt{aa}$: $\delta_{\rm CKM}$",      10),
        ("abbAB", r"$\mathtt{abbAB}$: $\theta_{12}^\nu$",   -12),
        ("B",     r"$\mathtt{b}$: $\approx\pi/2$",           10),
        ("aaabb", r"$\mathtt{aaabb}$: $\theta_{23}^{\rm CKM}$", -12),
    ],
}

# ── Layout: 2x2 grid (top=high angles, bottom=low angles) x 2 manifolds ──────
fig = plt.figure(figsize=(10, 7))

# GridSpec: top row taller (30-90 deg), bottom row shorter (0-20 deg)
from matplotlib.gridspec import GridSpec
gs = GridSpec(2, 2, height_ratios=[2.8, 1.2], hspace=0.06, wspace=0.08)

axes_top = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
axes_bot = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

DATASETS = [
    (m003, "m003", r"$M_{\rm PMNS}=\mathtt{m003}$", "#2ca02c"),
    (m006, "m006", r"$M_{\rm CKM}=\mathtt{m006}$",  "#ff7f0e"),
]

for col, (df, mname, title, color) in enumerate(DATASETS):
    ax_top = axes_top[col]
    ax_bot = axes_bot[col]

    rng = np.random.default_rng(42)
    jitter = rng.uniform(-0.10, 0.10, len(df))
    x = df["word_len"] + jitter
    y = df["phi_fold"]

    # ── top panel: 20–92 deg, linear ──────────────────────────────────────
    mask_top = y >= 20
    ax_top.scatter(x[mask_top], y[mask_top], color=color, s=26,
                   zorder=3, edgecolors="white", linewidths=0.3, alpha=0.92)
    # also plot small-angle points faintly to show they exist
    ax_top.scatter(x[~mask_top], y[~mask_top], color=color, s=8,
                   zorder=2, edgecolors="none", alpha=0.18)

    ax_top.set_ylim(18, 92)
    ax_top.set_yscale("linear")
    ax_top.set_xlim(0.5, 5.7)
    ax_top.set_xticks([1,2,3,4,5])
    ax_top.tick_params(bottom=False, labelbottom=False, labelsize=8)
    ax_top.grid(True, alpha=0.18, linewidth=0.5)
    ax_top.set_title(title, fontsize=11, pad=5)

    # SM lines in top panel
    for tval, tlabel, tc in SM:
        if tval >= 18:
            ax_top.axhline(tval, color=tc, alpha=0.35, lw=1.0, ls="--", zorder=1)
            ax_top.text(5.55, tval, tlabel, va="center", ha="left",
                        fontsize=7, color=tc, alpha=0.85)

    # ── bottom panel: 0–20 deg, log scale ─────────────────────────────────
    mask_bot = y <= 20
    # add small epsilon to avoid log(0)
    y_bot = y[mask_bot].clip(lower=0.15)
    ax_bot.scatter(x[mask_bot], y_bot, color=color, s=26,
                   zorder=3, edgecolors="white", linewidths=0.3, alpha=0.92)

    ax_bot.set_ylim(0.12, 20)
    ax_bot.set_yscale("log")
    ax_bot.set_xlim(0.5, 5.7)
    ax_bot.set_xticks([1,2,3,4,5])
    ax_bot.set_xlabel("Word length $|w|$", fontsize=10)
    ax_bot.tick_params(labelsize=8)
    ax_bot.grid(True, alpha=0.18, linewidth=0.5, which="both")

    # SM lines in bottom panel
    for tval, tlabel, tc in SM:
        if tval <= 20:
            ax_bot.axhline(tval, color=tc, alpha=0.35, lw=1.0, ls="--", zorder=1)
            ax_bot.text(5.55, tval, tlabel, va="center", ha="left",
                        fontsize=7, color=tc, alpha=0.85)

    # ── broken axis cosmetics ──────────────────────────────────────────────
    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    # diagonal break marks
    d = 0.012
    kwargs = dict(transform=ax_top.transAxes, color="k", clip_on=False, lw=0.8)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    ax_top.plot((1-d, 1+d), (-d, +d), **kwargs)
    kwargs.update(transform=ax_bot.transAxes)
    ax_bot.plot((-d, +d), (1-d, 1+d), **kwargs)
    ax_bot.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

    # ── annotations: top 4 hits ────────────────────────────────────────────
    for word_key, label, dy in ANNOT.get(mname, []):
        row = df[df["word"].str.upper() == word_key.upper()]
        if row.empty:
            row = df[df["word"] == word_key]
        if row.empty:
            continue
        xp = float(row["word_len"].values[0])
        yp = float(row["phi_fold"].values[0])
        target_ax = ax_top if yp >= 20 else ax_bot
        target_ax.annotate(
            label, xy=(xp, yp),
            xytext=(xp + 0.3, yp * (1 + dy/100) if target_ax is ax_bot else yp + dy*0.4),
            fontsize=6.8, color=color,
            arrowprops=dict(arrowstyle="-", color=color, lw=0.6, alpha=0.75),
            va="center",
        )

# ── shared y labels ───────────────────────────────────────────────────────────
fig.text(0.005, 0.72, r"$|\phi(\gamma)|$ (deg), linear",
         va="center", rotation="vertical", fontsize=9)
fig.text(0.005, 0.22, r"$|\phi(\gamma)|$ (deg), log",
         va="center", rotation="vertical", fontsize=9)

# ── legend ────────────────────────────────────────────────────────────────────
ckm_l  = mlines.Line2D([], [], color="#1f77b4", ls="--", alpha=0.7, label="CKM SM target")
pmns_l = mlines.Line2D([], [], color="#d62728", ls="--", alpha=0.7, label="PMNS SM target")
fig.legend(handles=[ckm_l, pmns_l], loc="lower center", ncol=2,
           fontsize=9, framealpha=0.85, bbox_to_anchor=(0.5, 0.0))

fig.suptitle(
    r"A-factor twist spectrum $|\phi(\gamma)|$ vs.\ word length "
    r"for $\mathtt{m003}$ and $\mathtt{m006}$",
    fontsize=11, y=1.005)

out = r"C:\dev\framework\papers\hyperbolic-flavor-twist\fig_twist_spectrum.pdf"
fig.savefig(out, format="pdf", bbox_inches="tight", dpi=300)
print(f"Saved: {out}")
