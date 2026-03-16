"""
paper8_analysis.py
Reads twist_census_m003.csv and twist_census_m006.csv,
deduplicates by phi (rounded to 3 decimal places),
and produces a ranked coincidence table for Paper 8.
"""

import pandas as pd
import numpy as np

SM_ANGLES = {
    "theta12_CKM":   13.04,
    "theta13_CKM":    0.201,
    "theta23_CKM":    2.38,
    "delta_CKM":     68.0,
    "theta12_NU":    33.41,
    "theta13_NU":     8.54,
    "theta23_NU":    49.1,
    "delta_CP":     197.0,
    "90deg":         90.0,
    "7xtheta12_CKM": 7 * 13.04,
}

SM_RATIOS = {
    "mb_mc":  3.2913,
    "MZ_MW":  1.1345,
}

ANGLE_TOL_DEG = 3.0
RATIO_TOL_PCT = 5.0

def fold_candidates(phi_deg):
    p = abs(phi_deg) % 360
    candidates = {p}
    if p > 180:
        candidates.add(360 - p)
    candidates.add(180 - p if p <= 180 else p - 180)
    return {c for c in candidates if 0 <= c <= 180}

def load_census(path, manifold_name):
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    df["manifold"] = manifold_name
    return df

m003 = load_census(r"C:\dev\hyperbolic-flavor-scan\twist_census_m003.csv", "m003")
m006 = load_census(r"C:\dev\hyperbolic-flavor-scan\twist_census_m006.csv", "m006")
df_all = pd.concat([m003, m006], ignore_index=True)
print(f"Loaded {len(m003)} m003 rows, {len(m006)} m006 rows")

df_all["phi_key"] = df_all["phi_deg"].round(3)
df_all["word_len"] = df_all["word"].str.len()
df_all = df_all.sort_values(["manifold", "phi_key", "word_len", "word"])
df_dedup = df_all.drop_duplicates(subset=["manifold", "phi_key"], keep="first").copy()
df_dedup = df_dedup.sort_values(["manifold", "phi_key"]).reset_index(drop=True)

print(f"\nAfter deduplication: {len(df_dedup)} unique (manifold, phi) pairs")
print(f"  m003: {(df_dedup.manifold == 'm003').sum()}")
print(f"  m006: {(df_dedup.manifold == 'm006').sum()}")

angle_rows = []
for _, row in df_dedup.iterrows():
    phi = row["phi_deg"]
    for cand in fold_candidates(phi):
        for tgt_name, tgt_val in SM_ANGLES.items():
            err = abs(cand - tgt_val)
            if err <= ANGLE_TOL_DEG:
                folding = "direct" if abs(phi - cand) < 0.01 else "180-phi"
                angle_rows.append({
                    "manifold":    row["manifold"],
                    "word":        row["word"].strip(),
                    "length":      int(row["word_len"]),
                    "phi_deg":     round(phi, 3),
                    "phi_used":    round(cand, 3),
                    "folding":     folding,
                    "target":      tgt_name,
                    "target_val":  tgt_val,
                    "abs_err":     round(err, 3),
                    "rel_err_pct": round(100 * err / tgt_val if tgt_val != 0 else 0, 2),
                    "type":        "angle",
                })

df_angles = pd.DataFrame(angle_rows)

ratio_rows = []
for manifold in ["m003", "m006"]:
    sub = df_dedup[df_dedup.manifold == manifold].sort_values("word_len")
    rows = sub.to_dict("records")
    for i, r1 in enumerate(rows):
        for r2 in rows[i+1:]:
            ml1 = r1.get("mod_lambda", np.nan)
            ml2 = r2.get("mod_lambda", np.nan)
            if np.isnan(ml1) or np.isnan(ml2) or ml2 == 0:
                continue
            for ratio, w1, w2 in [(ml1/ml2, r1["word"].strip(), r2["word"].strip()),
                                   (ml2/ml1, r2["word"].strip(), r1["word"].strip())]:
                if ratio < 1:
                    continue
                for tgt_name, tgt_val in SM_RATIOS.items():
                    err_pct = 100 * abs(ratio - tgt_val) / tgt_val
                    if err_pct <= RATIO_TOL_PCT:
                        ratio_rows.append({
                            "manifold":    manifold,
                            "word":        f"{w1}/{w2}",
                            "length":      max(len(w1), len(w2)),
                            "phi_deg":     None,
                            "phi_used":    None,
                            "folding":     "ratio",
                            "target":      tgt_name,
                            "target_val":  tgt_val,
                            "abs_err":     round(abs(ratio - tgt_val), 4),
                            "rel_err_pct": round(err_pct, 2),
                            "type":        "ratio",
                            "ratio_val":   round(ratio, 4),
                        })

df_ratios = pd.DataFrame(ratio_rows) if ratio_rows else pd.DataFrame()
df_combined = pd.concat([df_angles, df_ratios], ignore_index=True, sort=False)
df_combined = df_combined.sort_values("rel_err_pct").reset_index(drop=True)

headline_rows = []
for tgt in list(SM_ANGLES.keys()) + list(SM_RATIOS.keys()):
    sub = df_combined[df_combined.target == tgt]
    if sub.empty:
        continue
    headline_rows.append(sub.nsmallest(1, "rel_err_pct").iloc[0])

df_headline = pd.DataFrame(headline_rows).sort_values("rel_err_pct")

out_all  = r"C:\dev\hyperbolic-flavor-scan\paper8_coincidences.csv"
out_head = r"C:\dev\hyperbolic-flavor-scan\paper8_headline.csv"
df_combined.to_csv(out_all, index=False)
df_headline.to_csv(out_head, index=False)
print(f"\nSaved:\n  {out_all}\n  {out_head}")

print("\n" + "="*80)
print("HEADLINE COINCIDENCES")
print("="*80)
cols = ["manifold","word","length","phi_used","target","target_val","abs_err","rel_err_pct"]
print(df_headline[[c for c in cols if c in df_headline.columns]].to_string(index=False))

print("\n" + "="*80)
print("ALL SUB-1% MATCHES")
print("="*80)
sub1 = df_combined[df_combined.rel_err_pct < 1.0]
print(sub1[[c for c in cols if c in sub1.columns]].to_string(index=False))

print("\nDone.")
