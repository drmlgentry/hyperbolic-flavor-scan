import pandas as pd
import numpy as np

df6 = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len7_m006.csv")

# Deduplicate by phi_fold rounded to 4 dp
df6["phi_key"] = df6["phi_fold"].round(4)
df6 = df6.sort_values(["phi_key","length","word"])
dedup = df6.drop_duplicates("phi_key", keep="first").sort_values("phi_key")

print("Full deduplicated phi_fold spectrum on m006 (lengths 1-7):")
print(f"{'phi_fold':>12}  {'length':>6}  {'word':<12}  notes")
print("-"*55)
for _, row in dedup.iterrows():
    phi = row["phi_fold"]
    notes = ""
    if phi < 0.01:
        notes = "<-- near-identity"
    elif abs(phi - 0.201) < 0.05:
        notes = "<-- theta13_CKM TARGET"
    elif abs(phi - 2.38) < 0.3:
        notes = "<-- theta23_CKM"
    elif abs(phi - 68.0) < 1.0:
        notes = "<-- delta_CKM"
    elif abs(phi - 33.41) < 1.0:
        notes = "<-- theta12_nu"
    print(f"{phi:>12.6f}  {int(row['length']):>6}  {row['word']:<12}  {notes}")

print(f"\nSpectral gap analysis:")
vals = sorted(dedup["phi_fold"].values)
print(f"Smallest:  {vals[0]:.6f} deg")
print(f"Second:    {vals[1]:.6f} deg")
gaps = [(vals[i+1]-vals[i], vals[i], vals[i+1]) for i in range(len(vals)-1)]
gaps.sort(reverse=True)
print(f"\nLargest gaps:")
for gap, lo, hi in gaps[:5]:
    print(f"  [{lo:.4f} -> {hi:.4f}]  gap = {gap:.4f} deg")
print(f"\ntheta13_CKM = 0.201 deg falls in gap: "
      f"{[f'[{lo:.4f},{hi:.4f}]' for g,lo,hi in gaps if lo < 0.201 < hi]}")
