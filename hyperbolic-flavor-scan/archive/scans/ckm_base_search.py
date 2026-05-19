import snappy
import numpy as np

M2 = snappy.OrientableClosedCensus[43]
v2 = float(M2.volume())

print(f"M2/CKM: vol={v2:.8f}, H1={M2.homology()}")
print(f"Searching for base manifolds (checking first 500 census entries)...\n")

# For each degree, what volume would the base have?
for deg in [2, 3, 4, 5]:
    target = v2 / deg
    print(f"Degree-{deg} base would have volume: {target:.8f}")
    for i in range(500):
        M = snappy.OrientableClosedCensus[i]
        v = float(M.volume())
        if abs(v - target) < 1e-4:
            print(f"  MATCH: census[{i}] {M.name()} vol={v:.8f} H1={M.homology()}")

# Also: does anything in the census have volume = v2 * k for k=2,3?
print("\nSearching for manifolds that M2 covers (M2 as base):")
for i in range(500):
    M = snappy.OrientableClosedCensus[i]
    v = float(M.volume())
    for deg in [2, 3]:
        if abs(v - v2 * deg) < 1e-4:
            print(f"  census[{i}] {M.name()} vol={v:.8f} H1={M.homology()} -- degree-{deg} cover of M2")
