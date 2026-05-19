import snappy
import numpy as np
from math import gcd

M3 = snappy.OrientableClosedCensus[1]   # PMNS base
v3 = float(M3.volume())

print(f"M3/PMNS base: {M3.name()}, vol={v3:.8f}, H1={M3.homology()}")
print(f"Searching full OrientableClosedCensus ({len(snappy.OrientableClosedCensus)} entries)")
print(f"Looking for degree-2, 3, and 4 covers of M3\n")

deg2, deg3, deg4 = [], [], []

n = len(snappy.OrientableClosedCensus)
for i in range(n):
    if i % 2000 == 0:
        print(f"  Progress: {i}/{n}...")
    M = snappy.OrientableClosedCensus[i]
    v = float(M.volume())
    if abs(v - 2*v3) < 1e-4:
        deg2.append((i, M.name(), v, str(M.homology())))
    if abs(v - 3*v3) < 1e-4:
        deg3.append((i, M.name(), v, str(M.homology())))
    if abs(v - 4*v3) < 1e-4:
        deg4.append((i, M.name(), v, str(M.homology())))

print(f"\n=== Degree-2 covers of M3/PMNS (vol={2*v3:.6f}) ===")
print(f"Count: {len(deg2)}")
for idx, name, vol, h1 in deg2:
    print(f"  census[{idx:5d}] {name:10s} vol={vol:.8f} H1={h1}")

print(f"\n=== Degree-3 covers of M3/PMNS (vol={3*v3:.6f}) ===")
print(f"Count: {len(deg3)}")
for idx, name, vol, h1 in deg3:
    print(f"  census[{idx:5d}] {name:10s} vol={vol:.8f} H1={h1}")

print(f"\n=== Degree-4 covers of M3/PMNS (vol={4*v3:.6f}) ===")
print(f"Count: {len(deg4)}")
for idx, name, vol, h1 in deg4:
    print(f"  census[{idx:5d}] {name:10s} vol={vol:.8f} H1={h1}")

# Arithmetic predictions from filling slopes (-2,3) and (-5,2)
print(f"\n=== Arithmetic Predictions from Filling Slopes ===")
pmns = np.array([-2, 3])
ckm  = np.array([-5, 2])

print(f"PMNS slope: {pmns},  CKM slope: {ckm}")
print(f"Intersection (det): {abs(int(pmns[0]*ckm[1] - pmns[1]*ckm[0]))}  -> predicts prime 11")
print(f"Sum |p|: {abs(int(pmns[0]))+abs(int(ckm[0]))}                     -> predicts prime 7")
print(f"Sum |q|: {abs(int(pmns[1]))+abs(int(ckm[1]))}                     -> predicts prime 5 (base torsion)")

# What new invariants arise at degree 3?
# Degree-3 subgroups of pi_1 relate to 3*(-2,3) = (-6,9) ~ (-2,3) scaled
# or combinations. Natural candidates:
p, q = -2, 3
candidates_deg3 = [
    ("3|p|",          3*abs(p)),
    ("3|q|",          3*abs(q)),
    ("|p|+2|q|",      abs(p)+2*abs(q)),
    ("2|p|+|q|",      2*abs(p)+abs(q)),
    ("|p|*|q|",       abs(p)*abs(q)),
    ("|p|^2+|q|",     abs(p)**2+abs(q)),
    ("|q|^2-|p|",     abs(q)**2-abs(p)),
    ("|p|^2+|q|^2",   abs(p)**2+abs(q)**2),
    ("det*3",         3*11),
    ("det+7",         11+7),
]
print(f"\nPredicted new primes at degree 3 (from PMNS slope (-2,3)):")
for label, val in candidates_deg3:
    # check if prime
    is_prime = val > 1 and all(val % d != 0 for d in range(2, int(val**0.5)+1))
    flag = " <- PRIME" if is_prime else ""
    print(f"  {label:20s} = {val}{flag}")
