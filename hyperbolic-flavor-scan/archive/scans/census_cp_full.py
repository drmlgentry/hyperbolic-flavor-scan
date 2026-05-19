import snappy
import numpy as np

def get_h1_info(M):
    """Get torsion order and abelianization coefficients robustly."""
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
    
    # Build relation matrix
    R = np.zeros((len(relators), n), dtype=int)
    for ri, rel in enumerate(relators):
        for c in rel:
            g = c.lower()
            if g not in gens: continue
            idx = gens.index(g)
            R[ri, idx] += 1 if c.islower() else -1
    
    # Find abelianization: solve R*x = 0 mod p to find generator coefficients
    # Strategy: row reduce over Z, find column that achieves order p
    coeffs = np.zeros(n, dtype=int)
    found = False
    for ri in range(len(R)):
        for ci in range(n):
            val = R[ri, ci] % p
            if val == 0:
                continue
            # This row gives: val * gen_ci + (others) = 0 mod p
            # If val and p are coprime, gen_ci is determined by others
            try:
                inv_val = pow(int(val), -1, p)
                # Set coefficients: express all gens in terms of gen_0
                coeffs[ci] = 1
                for cj in range(n):
                    if cj != ci:
                        coeffs[cj] = (-R[ri, cj] * inv_val) % p
                found = True
                break
            except:
                continue
        if found:
            break
    
    if not found or all(c == 0 for c in coeffs):
        # Fallback: use elementary divisors structure
        # For Z/p with 2 generators, try simple assignment
        coeffs = np.array([1, 0], dtype=int)[:n]
    
    return coeffs, p, gens

def word_class(word, coeffs, gens, p):
    val = 0
    for c in word:
        g = c.lower()
        if g not in gens: continue
        idx = gens.index(g)
        val += int(coeffs[idx]) if c.islower() else -int(coeffs[idx])
    return val % p

def j_geom(k1, k2, k3, p):
    if k1==k2 or k2==k3 or k1==k3:
        return 0.0
    phi = [2*np.pi*k/p for k in [k1,k2,k3]]
    return abs(np.sin(phi[0]-phi[1]) +
               np.sin(phi[1]-phi[2]) +
               np.sin(phi[2]-phi[0]))

def max_j_for_p(p):
    """Theoretical maximum of j_geom for Z/p."""
    best = 0.0
    for k1 in range(p):
        for k2 in range(k1+1, p):
            for k3 in range(k2+1, p):
                j = j_geom(k1, k2, k3, p)
                if j > best:
                    best = j
    return best

# Build word list from 2 generators (generic)
def make_words(gens, max_len=3):
    words = []
    alpha = gens[0] if len(gens) > 0 else 'a'
    beta  = gens[1] if len(gens) > 1 else 'b'
    ALPHA, BETA = alpha.upper(), beta.upper()
    for l in range(1, max_len+1):
        import itertools
        for combo in itertools.product([alpha,ALPHA,beta,BETA], repeat=l):
            w = ''.join(combo)
            # Skip trivial cancellations
            if any(w[i].lower()==w[i+1].lower() and w[i]!=w[i+1]
                   for i in range(len(w)-1)):
                continue
            words.append(w)
    return words[:60]  # limit

print("Scanning census[1:500] for cyclic torsion H1=Z/p...")
print(f"{'Idx':>5} {'p':>4} {'n_gens':>6} {'J_max':>8} {'J/J_max':>8} {'Label':>14} {'H1':>8}")
print("-"*60)

all_results = []
for idx in range(1, 501):
    try:
        M = snappy.OrientableClosedCensus[idx]
        coeffs, p, gens = get_h1_info(M)
        if coeffs is None:
            continue
        
        words = make_words(gens, max_len=3)
        classes = {}
        for w in words:
            c = word_class(w, coeffs, gens, p)
            classes[w] = c
        
        # Find max J over all word triples
        wlist = list(classes.keys())
        j_max = 0.0
        best_triple = None
        for i in range(min(len(wlist),40)):
            for j in range(i+1, min(len(wlist),40)):
                for k in range(j+1, min(len(wlist),40)):
                    J = j_geom(classes[wlist[i]],
                               classes[wlist[j]],
                               classes[wlist[k]], p)
                    if J > j_max:
                        j_max = J
                        best_triple = (wlist[i],wlist[j],wlist[k],
                                       classes[wlist[i]],classes[wlist[j]],
                                       classes[wlist[k]])
        
        j_theory = max_j_for_p(p) if p <= 20 else max_j_for_p(5)
        ratio = j_max / j_theory if j_theory > 0 else 0
        
        label = ("DEGENERATE" if j_max < 0.01 else
                 "HIGH-CP" if ratio > 0.9 else "MED-CP")
        
        all_results.append((idx, p, len(gens), j_max, ratio,
                            label, str(M.homology()), best_triple))
        
    except:
        continue

print(f"\nTotal with cyclic Z/p torsion: {len(all_results)}")
from collections import Counter
labels = Counter(r[5] for r in all_results)
for k,v in sorted(labels.items()):
    print(f"  {k}: {v}")

print(f"\nBy torsion order:")
by_p = Counter(r[1] for r in all_results)
for p,cnt in sorted(by_p.items()):
    sub = [r for r in all_results if r[1]==p]
    j_vals = [r[3] for r in sub]
    print(f"  Z/{p:3d}: {cnt:4d} manifolds, "
          f"J_max range [{min(j_vals):.3f}, {max(j_vals):.3f}]")

print(f"\nZ/5 manifolds (our flavor manifolds):")
for r in [x for x in all_results if x[1]==5]:
    print(f"  idx={r[0]:4d} J_max={r[3]:.4f} ratio={r[4]:.3f} "
          f"best_triple={r[7][:3] if r[7] else 'none'} "
          f"classes={r[7][3:] if r[7] else 'none'}")
