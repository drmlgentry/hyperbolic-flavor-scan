import snappy
import numpy as np

# The correct approach: compute actual J from the geometric construction
# for all Z/5 manifolds, using the same word-triple selection method
# as the CKM paper

def get_coeffs_stable(M):
    """Stable coefficient computation using known relator structure."""
    h1 = M.homology()
    ed = h1.elementary_divisors()
    if len(ed) != 1:
        return None, None, None
    p = int(ed[0])
    if p < 3:
        return None, None, None
    
    fg = M.fundamental_group()
    gens = fg.generators()
    n = len(gens)
    relators = fg.relators()
    
    # Build relation matrix carefully
    R = []
    for rel in relators:
        row = [0]*n
        for c in rel:
            g = c.lower()
            if g not in gens: continue
            idx = gens.index(g)
            row[idx] += 1 if c.islower() else -1
        R.append(row)
    R = np.array(R, dtype=int)
    
    # Find the column that has a coefficient divisible by p
    # (this is the generator with torsion order p)
    for ri in range(len(R)):
        for ci in range(n):
            v = R[ri, ci] % p
            if v == 0 and R[ri, ci] != 0:
                # This row says p | R[ri,ci], so generator ci
                # has its order divided. Look at other row.
                pass
    
    # Canonical approach: use the Smith-like structure
    # For H1=Z/p with 2 generators a,b:
    # There must be a relator with coefficient p on one generator
    # That generator is the torsion generator (order p)
    # The other generator maps to some multiple k mod p
    
    torsion_col = None
    torsion_row = None
    for ri in range(len(R)):
        for ci in range(n):
            if abs(R[ri, ci]) == p:
                torsion_col = ci
                torsion_row = ri
                break
        if torsion_col is not None:
            break
    
    if torsion_col is None:
        return None, None, None
    
    # The OTHER generator maps to some element in Z/p
    # From another relator: sum of coefficients = 0 mod p
    # Use a relator that's different from torsion_row
    coeffs = [0]*n
    coeffs[torsion_col] = 1
    
    for ri in range(len(R)):
        if ri == torsion_row:
            continue
        # This relator: sum_j R[ri,j]*c[j] = 0 mod p
        # We know c[torsion_col] = 1
        # Solve for others
        rhs = -R[ri, torsion_col] % p
        for cj in range(n):
            if cj == torsion_col:
                continue
            if R[ri, cj] % p != 0:
                try:
                    inv = pow(int(R[ri, cj] % p), -1, p)
                    coeffs[cj] = (rhs * inv) % p
                    break
                except:
                    pass
        if any(coeffs[cj] != 0 for cj in range(n) if cj != torsion_col):
            break
    
    return coeffs, p, gens

def word_class(word, coeffs, gens, p):
    val = 0
    for c in word:
        g = c.lower()
        if g not in gens: continue
        idx = gens.index(g)
        val += int(coeffs[idx]) if c.islower() else -int(coeffs[idx])
    return val % p

# Verify on m003 and m006
print("=== Verifying coefficient solver ===")
for idx, name, expected in [
    (43, "m006 (CKM)", {"aaB":3, "AbA":2, "AAb":2}),
    (1,  "m003 (PMNS)", None)
]:
    M = snappy.OrientableClosedCensus[idx]
    coeffs, p, gens = get_coeffs_stable(M)
    if coeffs is None:
        print(f"{name}: FAILED")
        continue
    
    print(f"\n{name}: coeffs={coeffs}, p={p}, gens={gens}")
    
    fg = M.fundamental_group()
    print(f"  Relators: {fg.relators()}")
    
    for w in ["a","b","A","B","aaB","AbA","AAb","aab","abA","aBa"]:
        c = word_class(w, coeffs, gens, p)
        print(f"  [{w}] = {c}")
    
    if expected:
        ok = all(word_class(w,coeffs,gens,p)==v 
                 for w,v in expected.items())
        print(f"  Matches expected {expected}: {ok}")
