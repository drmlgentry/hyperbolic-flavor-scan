import snappy
import numpy as np

# Established: M3 (census[1]) -> M1 (census[39]) is degree-2 cover
# Questions:
# 1. Is there a degree-2 cover of M1 in the census?
# 2. Is M3 itself a cover of something smaller?
# 3. Does M2/CKM fit anywhere in a covering tower?

M2 = snappy.OrientableClosedCensus[43]
M3 = snappy.OrientableClosedCensus[1]

v2 = float(M2.volume())
v3 = float(M3.volume())

print(f"M2/CKM volume: {v2:.8f}")
print(f"M3/PMNS volume: {v3:.8f}")

# Search census for manifolds whose volume is v3/2, v3/3 (M3 as cover)
# or v2/2, v2/3 (M2 as cover)
print("\nSearching OrientableClosedCensus for volume relationships...")
print("(Checking first 200 manifolds)")

candidates = []
for i in range(200):
    M = snappy.OrientableClosedCensus[i]
    v = float(M.volume())
    
    # Is this a base that M3 covers?
    for deg in [2, 3, 4]:
        if abs(v * deg - v3) < 1e-4:
            candidates.append((i, M.name(), v, 'M3-covers-this', deg))
    
    # Is this a base that M2 covers?
    for deg in [2, 3, 4]:
        if abs(v * deg - v2) < 1e-4:
            candidates.append((i, M.name(), v, 'M2-covers-this', deg))
    
    # Is this a degree-2 cover of M3?
    if abs(v - 2*v3) < 1e-4 and i != 39:
        candidates.append((i, M.name(), v, 'covers-M3-deg2', 2))
    
    # Is this a degree-2 cover of M2?
    if abs(v - 2*v2) < 1e-4:
        candidates.append((i, M.name(), v, 'covers-M2-deg2', 2))

print("\nCandidates found:")
if candidates:
    for idx, name, vol, rel, deg in candidates:
        print(f"  census[{idx:3d}] {name:8s} vol={vol:.6f}  {rel}  degree={deg}")
else:
    print("  None in first 200 manifolds")

# Also check: what is the volume ratio M2/M3?
print(f"\nM2/M3 volume ratio: {v2/v3:.8f}")
print(f"M2/M3 ratio - 2: {v2/v3 - 2:.8f}")
print(f"Weeks manifold volume: 0.94270736")
print(f"M3 / Weeks ratio: {v3/0.94270736:.8f}")
