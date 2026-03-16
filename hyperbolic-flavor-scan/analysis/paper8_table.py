"""
paper8_table.py
Reads paper8_coincidences.csv and outputs:
  1. A LaTeX longtable block for the full coincidence table (Table 2)
  2. A corrected headline table (Table 1) with canonical words

Filters: angles with rel_err_pct < 5.0, ratios with rel_err_pct < 2.0
One row per (manifold, target, type) keeping best match only for Table 1,
all passing matches for Table 2.
"""

import pandas as pd
import numpy as np

df = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\paper8_coincidences.csv")

# Fix canonical words: use shorter pair for ratios
# For ratios, word is "W1/W2" -- ensure W1 is longer (numerator = larger modulus)
def canonical_ratio(word):
    if "/" not in str(word):
        return word
    parts = word.split("/")
    if len(parts) != 2:
        return word
    w1, w2 = parts
    # Keep numerator as longer word; tiebreak lex
    if len(w1) < len(w2):
        return f"{w2}/{w1}"
    return word

df["word_canon"] = df["word"].apply(lambda w: canonical_ratio(str(w).strip()))

# ── Table 1: one best row per SM target ──────────────────────────────────────
TARGET_ORDER = [
    "theta12_CKM", "theta23_CKM", "delta_CKM",
    "theta12_NU", "theta13_NU", "theta23_NU", "delta_CP",
    "mb_mc", "MZ_MW", "7xtheta12_CKM",
]

PRETTY = {
    "theta12_CKM":   r"$\theta_{12}^{\rm CKM}=13.04°$",
    "theta23_CKM":   r"$\theta_{23}^{\rm CKM}=2.38°$",
    "delta_CKM":     r"$\delta_{\rm CKM}=68.0°$",
    "theta12_NU":    r"$\theta_{12}^{\nu}=33.41°$",
    "theta13_NU":    r"$\theta_{13}^{\nu}=8.54°$",
    "theta23_NU":    r"$\theta_{23}^{\nu}=49.1°$",
    "delta_CP":      r"$\delta_{\rm CP}=197°$",
    "mb_mc":         r"$m_b/m_c=3.291$",
    "MZ_MW":         r"$M_Z/M_W=1.1345$",
    "7xtheta12_CKM": r"$7\theta_{12}^{\rm CKM}=91.28°$",
}

FOLD_PRETTY = {
    "direct":  "direct",
    "180-phi": r"$180°-|\phi|$",
    "ratio":   "ratio",
}

print("="*70)
print("TABLE 1: HEADLINE COINCIDENCES (LaTeX rows)")
print("="*70)
print()

tab1_lines = []
for tgt in TARGET_ORDER:
    sub = df[df["target"] == tgt]
    if sub.empty:
        continue
    best = sub.nsmallest(1, "rel_err_pct").iloc[0]
    mfld   = best["manifold"]
    word   = best["word_canon"]
    length = int(best["length"]) if not pd.isna(best["length"]) else "?"
    phi    = f"{best['phi_used']:.3f}°" if not pd.isna(best.get("phi_used", float("nan"))) else "---"
    fold   = FOLD_PRETTY.get(str(best.get("folding","")).strip(), str(best.get("folding","")))
    tval   = best["target_val"]
    aerr   = best["abs_err"]
    rerr   = best["rel_err_pct"]
    label  = PRETTY.get(tgt, tgt)

    # Format errors
    if best.get("type") == "ratio":
        err_str = f"{aerr:.4f} ({rerr:.2f}\\%)"
        phi_str = "---"
        fold_str = "modulus ratio"
    else:
        err_str = f"{aerr:.3f}° ({rerr:.2f}\\%)"
        phi_str = phi
        fold_str = fold

    line = (f"    \\mfld{{{mfld}}} & \\word{{{word}}} & {length} & "
            f"{phi_str} & {fold_str} & {label} & {err_str} \\\\")
    tab1_lines.append(line)
    print(line)

# ── Table 2: full longtable ───────────────────────────────────────────────────
print()
print("="*70)
print("TABLE 2: FULL COINCIDENCES (LaTeX longtable body)")
print("="*70)
print()

# Filter
mask_angle = (df["type"] == "angle") & (df["rel_err_pct"] < 5.0)
mask_ratio = (df["type"] == "ratio") & (df["rel_err_pct"] < 2.0)
df_tab2 = df[mask_angle | mask_ratio].copy()
df_tab2 = df_tab2.sort_values(["target", "manifold", "rel_err_pct"])

prev_tgt = None
for _, row in df_tab2.iterrows():
    tgt = row["target"]
    if tgt != prev_tgt:
        print(f"    \\midrule")
        print(f"    \\multicolumn{{7}}{{l}}{{\\textit{{Target: {PRETTY.get(tgt, tgt)}}}}} \\\\")
        prev_tgt = tgt

    mfld  = row["manifold"]
    word  = row["word_canon"]
    leng  = int(row["length"]) if not pd.isna(row["length"]) else "?"
    rerr  = row["rel_err_pct"]

    if row["type"] == "ratio":
        phi_str  = "---"
        fold_str = "ratio"
    else:
        phi_val = row.get("phi_used", float("nan"))
        phi_str  = f"{phi_val:.3f}°" if not pd.isna(phi_val) else "---"
        fold_str = FOLD_PRETTY.get(str(row.get("folding","")).strip(), "")

    aerr = row["abs_err"]
    if row["type"] == "ratio":
        err_str = f"{aerr:.4f} ({rerr:.2f}\\%)"
    else:
        err_str = f"{aerr:.3f}° ({rerr:.2f}\\%)"

    print(f"    \\mfld{{{mfld}}} & \\word{{{word}}} & {leng} & "
          f"{phi_str} & {fold_str} & {err_str} \\\\")

# ── Save outputs ──────────────────────────────────────────────────────────────
out = r"C:\dev\hyperbolic-flavor-scan\paper8_tables.txt"
import sys
# Re-run capturing to file -- just print notice here
print()
print(f"Copy the output above into gentry-hyperbolic-flavor-twist.tex")
print(f"(Table 1 rows replace the tabular body; Table 2 rows go in the longtable stub)")
