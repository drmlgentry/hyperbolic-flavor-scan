import pandas as pd
import numpy as np

# PDG CKM moduli
V_CKM = np.array([
    [0.97401, 0.22650, 0.00361],
    [0.22650, 0.97359, 0.04153],
    [0.00854, 0.04053, 0.99910]
])

# Read the CSV file
df = pd.read_csv('scan_results.csv')

# Ensure fitness column is numeric
df['fitness'] = pd.to_numeric(df['fitness'], errors='coerce')
df = df.dropna(subset=['fitness'])

# Sort by fitness and take top 5
top5 = df.nsmallest(5, 'fitness')

print("Top 5 candidates (by fitness):")
for idx, row in top5.iterrows():
    print(f"\n{row['name']} (vol={row['volume']:.4f}, fitness={row['fitness']:.4f}, J={row['jarlskog']:.2e})")
    # Extract the 3x3 moduli matrix (columns u11..u33)
    U_mod = np.array([[row['u11'], row['u12'], row['u13']],
                      [row['u21'], row['u22'], row['u23']],
                      [row['u31'], row['u32'], row['u33']]])
    print("Predicted moduli matrix:")
    print(np.array_str(U_mod, precision=4, suppress_small=True))
    print("CKM target:")
    print(np.array_str(V_CKM, precision=4, suppress_small=True))
    diff = U_mod - V_CKM
    print("Difference:")
    print(np.array_str(diff, precision=4, suppress_small=True))