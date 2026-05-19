import snappy
import numpy as np

M2 = snappy.OrientableClosedCensus[43]
v2 = float(M2.volume())
target_deg2 = v2 / 2.0  # 1.01442655

print(f"M2/CKM vol={v2:.8f}, H1={M2.homology()}")
print(f"Searching full OrientableClosedCensus ({len(snappy.OrientableClosedCensus)} entries)")
print(f"Looking for degree-2 base at vol={target_deg2:.8f} +/- 1e-4\n")

matches = []
n = len(snappy.OrientableClosedCensus)
for i in range(n):
    if i % 1000 == 0:
        print(f"  Progress: {i}/{n}...")
    M = snappy.OrientableClosedCensus[i]
    v = float(M.volume())
    if abs(v - target_deg2) < 1e-4:
        matches.append((i, M.name(), v, str(M.homology())))

print(f"\nDegree-2 base candidates (vol~{target_deg2:.6f}):")
if matches:
    for idx, name, vol, h1 in matches:
        print(f"  census[{idx:5d}] {name:10s} vol={vol:.8f} H1={h1}")
else:
    print("  None found in full census.")

# Also check degree-3 base while we're here
target_deg3 = v2 / 3.0
print(f"\nDegree-3 base candidates (vol~{target_deg3:.6f}):")
matches3 = []
for i in range(n):
    M = snappy.OrientableClosedCensus[i]
    v = float(M.volume())
    if abs(v - target_deg3) < 1e-4:
        matches3.append((i, M.name(), v, str(M.homology())))
if matches3:
    for idx, name, vol, h1 in matches3:
        print(f"  census[{idx:5d}] {name:10s} vol={vol:.8f} H1={h1}")
else:
    print("  None found in full census.")

# Bonus: report M2's systole and chirality
ls2 = M2.length_spectrum(3)
systole = min([float(abs(x.length.real)) for x in ls2])
print(f"\nM2 systole (shortest geodesic): {systole:.8f}")
print(f"M2 is orientable: {M2.is_orientable()}")
try:
    print(f"M2 chirality (is amphicheiral): {M2.is_amphicheiral()}")
except:
    print("M2 chirality: not available via this method")
