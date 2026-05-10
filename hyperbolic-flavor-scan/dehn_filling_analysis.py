import snappy
import numpy as np

# The two cusped ancestors
m003 = snappy.Manifold("m003")
m006 = snappy.Manifold("m006")

print("=== Cusped Ancestors ===")
print(f"m003: vol={float(m003.volume()):.8f}, cusps={m003.num_cusps()}, H1={m003.homology()}")
print(f"m006: vol={float(m006.volume()):.8f}, cusps={m006.num_cusps()}, H1={m006.homology()}")

# Verify the Dehn fillings
print("\n=== Verifying Dehn Fillings ===")
m003_filled = m003.copy()
m003_filled.dehn_fill((-2, 3))
vol_filled = float(m003_filled.volume())
print(f"m003(-2,3): vol={vol_filled:.8f}  (M3/PMNS vol=0.98136883, diff={abs(vol_filled-0.98136883):.2e})")

m006_filled = m006.copy()
m006_filled.dehn_fill((-5, 2))
vol_filled2 = float(m006_filled.volume())
print(f"m006(-5,2): vol={vol_filled2:.8f}  (M2/CKM vol=2.02885309, diff={abs(vol_filled2-2.02885309):.2e})")

# Explore neighboring Dehn fillings of m003
# What other closed manifolds live near (-2,3)?
print("\n=== Neighboring Dehn Fillings of m003 (near PMNS point) ===")
print(f"{'(p,q)':12s}  {'volume':12s}  {'H1':20s}")
for dp in range(-2, 3):
    for dq in range(-2, 3):
        p = -2 + dp
        q =  3 + dq
        if abs(p) + abs(q) == 0:
            continue
        try:
            M = m003.copy()
            M.dehn_fill((p, q))
            v = float(M.volume())
            h = str(M.homology())
            print(f"  ({p:2d},{q:2d})      {v:.8f}    {h}")
        except:
            pass

# Explore neighboring Dehn fillings of m006
print("\n=== Neighboring Dehn Fillings of m006 (near CKM point) ===")
print(f"{'(p,q)':12s}  {'volume':12s}  {'H1':20s}")
for dp in range(-2, 3):
    for dq in range(-2, 3):
        p = -5 + dp
        q =  2 + dq
        if abs(p) + abs(q) == 0:
            continue
        try:
            M = m006.copy()
            M.dehn_fill((p, q))
            v = float(M.volume())
            h = str(M.homology())
            print(f"  ({p:2d},{q:2d})      {v:.8f}    {h}")
        except:
            pass

# The key question: is there any filling of m003 that gives M2/CKM volume?
print("\n=== Searching m003 fillings for CKM volume (2.02885309) ===")
target = 2.02885309
for p in range(-10, 11):
    for q in range(-10, 11):
        if p == 0 and q == 0:
            continue
        from math import gcd
        if abs(p) > 0 and abs(q) > 0 and gcd(abs(p), abs(q)) != 1:
            continue
        try:
            M = m003.copy()
            M.dehn_fill((p, q))
            v = float(M.volume())
            if abs(v - target) < 1e-4:
                print(f"  m003({p},{q}): vol={v:.8f}  H1={M.homology()}")
        except:
            pass

print("\n=== Searching m006 fillings for PMNS volume (0.98136883) ===")
target2 = 0.98136883
for p in range(-10, 11):
    for q in range(-10, 11):
        if p == 0 and q == 0:
            continue
        from math import gcd
        if abs(p) > 0 and abs(q) > 0 and gcd(abs(p), abs(q)) != 1:
            continue
        try:
            M = m006.copy()
            M.dehn_fill((p, q))
            v = float(M.volume())
            if abs(v - target2) < 1e-4:
                print(f"  m006({p},{q}): vol={v:.8f}  H1={M.homology()}")
        except:
            pass
