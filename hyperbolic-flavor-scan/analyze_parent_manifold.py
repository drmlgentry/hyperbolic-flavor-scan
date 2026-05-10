import pandas as pd
import numpy as np
import glob, os

# ── Load all scan CSVs ──────────────────────────────────────────────────────
csv_files = glob.glob("*.csv")
print(f"Found {len(csv_files)} CSV files: {csv_files}\n")

dfs = []
for f in csv_files:
    try:
        df = pd.read_csv(f)
        df["_source"] = f
        dfs.append(df)
        print(f"  {f}: {len(df)} rows, columns: {list(df.columns)}")
    except Exception as e:
        print(f"  {f}: FAILED ({e})")

if not dfs:
    print("No data loaded.")
    exit()

df = pd.concat(dfs, ignore_index=True)
print(f"\nTotal rows: {len(df)}")
print(f"Columns: {list(df.columns)}\n")

# ── Known anchor values ─────────────────────────────────────────────────────
VOL_m003 = 0.98136857679
VOL_m006 = 2.02988321282
VOL_sum  = VOL_m003 + VOL_m006          # 3.011
VOL_2m003 = 2 * VOL_m003               # 1.963 (idx39, known)
PHI = (1 + 5**0.5) / 2

print("=== TARGET VOLUMES TO SEARCH FOR ===")
print(f"  vol(m003)              = {VOL_m003:.6f}")
print(f"  vol(m006)              = {VOL_m006:.6f}")
print(f"  vol(idx39) = 2*m003    = {VOL_2m003:.6f}  [known degree-2 cover of m003]")
print(f"  vol(m003)+vol(m006)    = {VOL_sum:.6f}  [parent candidate]")
print(f"  phi*vol(m003)          = {PHI*VOL_m003:.6f}")
print(f"  phi*vol(m006)          = {PHI*VOL_m006:.6f}")
print(f"  2*vol(m006)            = {2*VOL_m006:.6f}\n")

# ── Detect volume column ─────────────────────────────────────────────────────
vol_col = None
for candidate in ["volume", "vol", "Volume", "manifold_volume"]:
    if candidate in df.columns:
        vol_col = candidate
        break
if vol_col is None:
    print("ERROR: No volume column found. Columns are:", list(df.columns))
    exit()
print(f"Using volume column: '{vol_col}'\n")

df[vol_col] = pd.to_numeric(df[vol_col], errors="coerce")
df = df.dropna(subset=[vol_col])

# ── Detect H1 column ─────────────────────────────────────────────────────────
h1_col = None
for candidate in ["H1", "h1", "homology", "first_homology", "H_1"]:
    if candidate in df.columns:
        h1_col = candidate
        break

# ── Search for parent manifold candidates ───────────────────────────────────
tol = 0.0001
targets = {
    "vol(m003)+vol(m006)": VOL_sum,
    "2*vol(m006)":         2 * VOL_m006,
    "phi*vol(m003)":       PHI * VOL_m003,
    "phi*vol(m006)":       PHI * VOL_m006,
    "3*vol(m003)":         3 * VOL_m003,
}

print("=== PARENT MANIFOLD SEARCH ===")
found_any = False
for label, target in targets.items():
    hits = df[np.abs(df[vol_col] - target) < tol]
    if len(hits) > 0:
        found_any = True
        print(f"\n  TARGET {label} = {target:.6f}  -->  {len(hits)} HIT(S):")
        cols_to_show = [vol_col]
        for c in ["manifold", "index", "census_index", "H1", "h1", "homology", "fitness", "words"]:
            if c in hits.columns:
                cols_to_show.append(c)
        print(hits[cols_to_show].to_string(index=False))
    else:
        print(f"  TARGET {label} = {target:.6f}  -->  no hits within tol={tol}")

if not found_any:
    print("\n  No exact parent candidates found. Showing 10 manifolds closest to vol(m003)+vol(m006):")
    df["_dist"] = np.abs(df[vol_col] - VOL_sum)
    near = df.nsmallest(10, "_dist")
    cols_to_show = [vol_col, "_dist"]
    for c in ["manifold", "index", "H1", "h1", "homology"]:
        if c in near.columns:
            cols_to_show.append(c)
    print(near[cols_to_show].to_string(index=False))

# ── H1=Z/5 manifolds already in scan ────────────────────────────────────────
if h1_col:
    print(f"\n=== ALL MANIFOLDS WITH H1 CONTAINING Z/5 (column: {h1_col}) ===")
    z5 = df[df[h1_col].astype(str).str.contains("5", na=False)]
    print(f"  {len(z5)} rows match. Unique H1 values: {z5[h1_col].unique()[:20]}")
    # Show those with H1 = Z/5 exactly
    exact = z5[z5[h1_col].astype(str).str.strip() == "Z/5"]
    print(f"  Exact Z/5: {len(exact)} manifolds")
    if len(exact) > 0:
        print(exact[[vol_col, h1_col] + (["manifold"] if "manifold" in exact.columns else [])].head(20).to_string(index=False))
else:
    print("\n  No H1 column found in data — cannot filter by torsion.")

# ── Volume distribution summary ──────────────────────────────────────────────
print(f"\n=== VOLUME RANGE IN SCAN DATA ===")
print(f"  min={df[vol_col].min():.4f}  max={df[vol_col].max():.4f}  "
      f"median={df[vol_col].median():.4f}  n={len(df)}")
print(f"  Fraction with vol < 2.0:  {(df[vol_col]<2.0).mean():.1%}")
print(f"  Fraction with vol < 3.1:  {(df[vol_col]<3.1).mean():.1%}")
print(f"  Fraction with vol < 4.0:  {(df[vol_col]<4.0).mean():.1%}")

print("\nDone.")
