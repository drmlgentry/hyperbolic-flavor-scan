import pandas as pd, numpy as np

df = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len10_m006.csv")
df["h1_class"] = df["word"].apply(lambda w: sum(1 if c.islower() else -1 for c in w) % 5)

# Filter out trivial/near-trivial PSL elements (|lambda| < 1.01 means lambda ~ +-1)
df_geo = df[df.mod_lambda > 1.01].copy()
n_removed = len(df) - len(df_geo)
print(f"Removed {n_removed} trivial/near-trivial elements (|lambda| <= 1.01)")
print(f"Remaining: {len(df_geo)} genuine loxodromic geodesics\n")

print("Spectral floors per class (genuine geodesics only, length <= 10):")
for k in range(5):
    sub = df_geo[df_geo.h1_class == k].sort_values("phi_fold")
    if len(sub):
        row = sub.iloc[0]
        print(f"  Class [{k}]: floor={row.phi_fold:.6f} deg  word={row.word}  len={int(row.length)}  |lam|={row.mod_lambda:.4f}")
    else:
        print(f"  Class [{k}]: empty")

# Gap analysis
floors = {}
for k in range(5):
    sub = df_geo[df_geo.h1_class == k]
    if len(sub): floors[k] = sub.phi_fold.min()

sorted_floors = sorted(floors.values())
print(f"\nGap analysis (genuine geodesics only):")
print(f"  All floors: {[f'{v:.4f}' for v in sorted_floors]}")
print(f"  Gap: [{sorted_floors[0]:.4f}, {sorted_floors[1]:.4f}] deg")
print(f"  theta13_CKM = 0.201 in gap? {sorted_floors[0] < 0.201 < sorted_floors[1]}")
print(f"  Max/min ratio: {sorted_floors[-1]/sorted_floors[0]:.0f}x")

print(f"\nFirst near-zero (phi_fold < 0.01) per class, genuine only:")
for k in range(5):
    sub = df_geo[(df_geo.h1_class==k) & (df_geo.phi_fold < 0.01)].sort_values("length")
    if len(sub):
        row = sub.iloc[0]
        print(f"  Class {k}: length {int(row.length)}, word={row.word}, phi_fold={row.phi_fold:.6f}, |lam|={row.mod_lambda:.4f}")
    else:
        print(f"  Class {k}: no genuine near-zero at length <= 10")
