import snappy
import numpy as np
from collections import Counter

def h1_class_map(M):
    """Get the map from generators to Z/p classes via relation matrix."""
    fg = M.fundamental_group()
    gens = fg.generators()
    relators = fg.relators()
    n = len(gens)
    
    # Build relation matrix
    R = []
    for rel in relators:
        row = [0]*n
        for c in rel:
            idx = gens.index(c.lower())
            row[idx] += 1 if c.islower() else -1
        R.append(row)
    R = np.array(R, dtype=int)
    
    # Find torsion order p from relation matrix
    h1 = M.homology()
    ed = h1.elementary_divisors()
    if len(ed) != 1:
        return None, None  # skip non-cyclic
    p = int(ed[0])
    if p < 3:
        return None, None
    
    # Solve: find coefficients (c_a, c_b) such that
    # a = c_a * generator, b = c_b * generator mod p
    # Use: from relator 0*a + p*b = 0 => b is generator (coeff 1)
    # and from relator 1*a + r*b = 0 => a = -r*b mod p
    coeffs = [0]*n
    # Find which generator has order p (diagonal entry = p)
    for i in range(n):
        if abs(R[0,i]) == p:
            # generator i has order p, coefficient = 1
            coeffs[i] = 1
            # solve for others
            for j in range(n):
                if j != i:
                    # from relator: R[1,i]*coeffs[i] + R[1,j]*coeffs[j] = 0 mod p
                    # coeffs[j] = -R[1,i]/R[1,j] mod p
                    if R[1,j] != 0:
                        inv_Rj = pow(int(R[1,j]), -1, p)
                        coeffs[j] = (-R[1,i] * inv_Rj) % p
            break
    
    return coeffs, p

def word_class(word, coeffs, p):
    """Compute homology class of word."""
    val = 0
    gens = 'ab'  # assume 2 generators
    for c in word:
        g = c.lower()
        idx = gens.index(g) if g in gens else -1
        if idx >= 0:
            val += coeffs[idx] if c.islower() else -coeffs[idx]
    return val % p

def j_geom(k1, k2, k3, p):
    """Geometric CP proxy."""
    if k1==k2 or k2==k3 or k1==k3:
        return 0.0
    phi = [2*np.pi*k/p for k in [k1,k2,k3]]
    return (np.sin(phi[0]-phi[1]) +
            np.sin(phi[1]-phi[2]) +
            np.sin(phi[2]-phi[0]))

print("Census CP classification (indices 1-300)")
print(f"{'Idx':>5} {'p':>4} {'H1':>8} {'Best J_geom':>12} {'Triple classes':>20} {'Type':>8}")
print("-"*65)

results = []
for idx in range(1, 301):
    try:
        M = snappy.OrientableClosedCensus[idx]
        coeffs, p = h1_class_map(M)
        if coeffs is None:
            continue
        
        # Try all triples of distinct classes mod p
        best_J = 0.0
        best_triple = None
        for k1 in range(p):
            for k2 in range(k1+1, p):
                for k3 in range(k2+1, p):
                    J = abs(j_geom(k1,k2,k3,p))
                    if J > best_J:
                        best_J = J
                        best_triple = (k1,k2,k3)
        
        # Check if optimal geodesic triple has degeneracy
        # For m006 (idx 43): we know {3,2,2} -> degenerate
        mtype = "CKM-like" if best_J < 0.1 else "PMNS-like"
        
        results.append((idx, p, str(M.homology()), best_J, best_triple, mtype))
        
    except Exception:
        continue

# Show interesting cases
results.sort(key=lambda x: -x[3])
print("Top CP-active manifolds:")
for r in results[:15]:
    print(f"{r[0]:>5} {r[1]:>4} {r[2]:>8} {r[3]:>12.4f} {str(r[4]):>20} {r[5]:>8}")

# Summary stats
ckm_like = [r for r in results if r[1]==5]
print(f"\nZ/5 manifolds: {len(ckm_like)}")
print(f"Max J_geom for Z/5: {max(r[3] for r in ckm_like):.4f}")
print(f"sin(2pi/5) = {np.sin(2*np.pi/5):.4f}")
