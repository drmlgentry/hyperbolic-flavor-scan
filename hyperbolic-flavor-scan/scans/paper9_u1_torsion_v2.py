"""
paper9_u1_torsion_v2.py
Refined analysis: does Z/5 torsion control spectral FLOOR per coset?
Key question: is the near-identity geodesic isolated in its coset?
"""
import snappy, numpy as np, pandas as pd

def homology_class(word, n=5):
    return sum(1 if c.islower() else -1 for c in word) % n

for idx, name in [(1,"m003"),(43,"m006")]:
    M = snappy.OrientableClosedCensus[idx]
    print(f"\n{'='*55}")
    print(f"{name}  H1={M.homology()}")

    df = pd.read_csv(f"C:\\dev\\hyperbolic-flavor-scan\\data\\twist_census_len7_{name}.csv")
    df["phi_fold"] = df["phi_deg"].apply(lambda p: min(abs(p)%180, 180-abs(p)%180))
    df["phi_key"]  = df["phi_fold"].round(3)
    df = df.sort_values(["phi_key","length","word"])
    df = df.drop_duplicates("phi_key", keep="first")
    df["h1_class"] = df["word"].apply(homology_class)

    print(f"\nSpectral floor per coset (min phi_fold):")
    floors = {}
    for coset in range(5):
        sub = df[df.h1_class == coset]
        floor = sub.phi_fold.min()
        floors[coset] = floor
        print(f"  Coset [{coset}]: min={floor:.4f} deg  "
              f"(word: {sub.loc[sub.phi_fold.idxmin(),'word']}  len={sub.loc[sub.phi_fold.idxmin(),'length']})")

    # Is there a systematic pattern?
    print(f"\nFloors sorted: {sorted(floors.values())}")
    print(f"Ratio max_floor/min_floor: {max(floors.values())/min(floors.values()):.1f}x")

    # For m006: does the near-identity belong to a specific coset?
    low = df[df.phi_fold < 1.0]
    if len(low):
        print(f"\nAll phi_fold < 1.0 deg:")
        print(low[["word","phi_fold","h1_class","length"]].to_string())
        print(f"All in same coset? {low.h1_class.nunique() == 1}")

    # Check: words in identity coset [0] -- what is their floor?
    id_coset = df[df.h1_class == 0]
    print(f"\nIdentity coset [0] floor: {id_coset.phi_fold.min():.4f} deg")
    print(f"  (word: {id_coset.loc[id_coset.phi_fold.idxmin(),'word']})")

    # 5x torsion pattern
    floor_vals = sorted(floors.values())
    print(f"\n5x pattern check:")
    for f in floor_vals:
        print(f"  5 x {f:.4f} = {5*f:.4f} deg")
