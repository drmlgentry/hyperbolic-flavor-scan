import snappy
import numpy as np
import pandas as pd

PHI = (1 + 5**0.5) / 2
VOL_m003 = 0.98136857679
VOL_m006 = 2.02988321282

targets = {
    "vol(m003)+vol(m006)": VOL_m003 + VOL_m006,
    "2*vol(m006)":         2 * VOL_m006,
    "3*vol(m003)":         3 * VOL_m003,
    "phi*vol(m006)":       PHI * VOL_m006,
    "phi*vol(m003)":       PHI * VOL_m003,
    "2*vol(m003)+vol(m006)": 2*VOL_m003 + VOL_m006,
}

print("Target volumes:")
for k,v in targets.items():
    print(f"  {k:35s} = {v:.8f}")
print()

# Scan OrientableClosedCensus
census = snappy.OrientableClosedCensus
n = len(census)
print(f"Census size: {n} manifolds\n")

rows = []
tol = 0.0005  # generous tolerance — volumes are exact to 6+ decimals

for i, M in enumerate(census):
    try:
        vol = M.volume()
        for label, target in targets.items():
            if abs(vol - target) < tol:
                try:
                    h1 = M.homology().__repr__()
                except:
                    h1 = "unknown"
                rows.append({
                    "census_index": i,
                    "name": M.name(),
                    "volume": vol,
                    "H1": h1,
                    "target_label": label,
                    "vol_diff": vol - target,
                })
                print(f"  HIT  idx={i:4d}  vol={vol:.8f}  H1={h1:20s}  target={label}")
    except Exception as e:
        pass

    if i % 500 == 0:
        print(f"  ... {i}/{n} scanned")

print(f"\nTotal hits: {len(rows)}")

df = pd.DataFrame(rows)
if len(df) > 0:
    out = "data/parent_manifold_search.csv"
    df.to_csv(out, index=False)
    print(f"Saved to {out}")
    print(df.to_string(index=False))
else:
    print("No hits found at tolerance 0.0005.")
    # Show the 10 closest to vol(m003)+vol(m006) regardless
    print("\nFinding 10 closest to vol(m003)+vol(m006) across full census...")
    target = VOL_m003 + VOL_m006
    near = []
    for i, M in enumerate(census):
        try:
            vol = M.volume()
            near.append((abs(vol-target), i, M.name(), vol))
        except:
            pass
    near.sort()
    print(f"{'rank':>4}  {'idx':>5}  {'name':>12}  {'volume':>12}  {'|diff|':>10}")
    for rank, (diff, idx, name, vol) in enumerate(near[:10]):
        print(f"{rank+1:4d}  {idx:5d}  {name:>12}  {vol:12.8f}  {diff:10.8f}")
