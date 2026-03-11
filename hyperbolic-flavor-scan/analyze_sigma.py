import pandas as pd
import numpy as np
import glob

sigmas = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
target_manifold = "m007"   # change to any name, e.g., "m011", "m004"

best_score = float('inf')
best_sigma = None
best_row = None

for s in sigmas:
    fname = f"scan_results_sigma_{s}.csv"
    try:
        df = pd.read_csv(fname)
        df['fitness'] = pd.to_numeric(df['fitness'], errors='coerce')
        df = df.dropna(subset=['fitness'])
        mask = df['name'] == target_manifold
        if mask.any():
            row = df.loc[mask].iloc[0]   # first occurrence (should be unique)
            score = row['fitness']
            print(f"sigma={s}: {target_manifold} fitness = {score:.4f}")
            if score < best_score:
                best_score = score
                best_sigma = s
                best_row = row
    except FileNotFoundError:
        print(f"sigma={s}: file not found")

if best_sigma is not None:
    print(f"\nBest sigma for {target_manifold}: {best_sigma} (fitness={best_score:.4f})")
    print("\nMixing matrix moduli at best sigma:")
    U_mod = np.array([[best_row['u11'], best_row['u12'], best_row['u13']],
                      [best_row['u21'], best_row['u22'], best_row['u23']],
                      [best_row['u31'], best_row['u32'], best_row['u33']]])
    print(np.array_str(U_mod, precision=4, suppress_small=True))
else:
    print(f"No data found for {target_manifold}")