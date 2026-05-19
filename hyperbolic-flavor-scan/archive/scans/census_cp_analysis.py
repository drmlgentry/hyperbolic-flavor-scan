import snappy
import numpy as np

# Deep analysis of the MED-CP cases and the Z/5 family
# Focus: what distinguishes manifolds with low J_max from high J_max?

def get_h1_info(M):
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
    R = np.zeros((len(relators), n), dtype=int)
    for ri, rel in enumerate(relators):
        for c in rel:
            g = c.lower()
            if g not in gens: continue
            idx = gens.index(g)
            R[ri, idx] += 1 if c.islower() else -1
    coeffs = np.zeros(n, dtype=int)
    for ri in range(len(R)):
        for ci in range(n):
            val = R[ri, ci] % p
            if val == 0: continue
            try:
                inv = pow(int(val), -1, p)
                coeffs[ci] = 1
                for cj in range(n):
                    if cj != ci:
                        coeffs[cj] = (-R[ri, cj] * inv) % p
                return coeffs, p, gens
            except:
                continue
    return np.array([1,0]), p, gens

def word_class(word, coeffs, gens, p):
    val = 0
    for c in word:
        g = c.lower()
        if g not in gens: continue
        val += int(coeffs[gens.index(g)]) if c.islower() else -int(coeffs[gens.index(g)])
    return val % p

def j_geom(k1,k2,k3,p):
    if k1==k2 or k2==k3 or k1==k3: return 0.0
    phi = [2*np.pi*k/p for k in [k1,k2,k3]]
    return abs(np.sin(phi[0]-phi[1])+np.sin(phi[1]-phi[2])+np.sin(phi[2]-phi[0]))

# For Z/5 manifolds: verify the CKM triple classes exactly
print("=== Z/5 manifolds: CKM triple (aaB, AbA, AAb) class analysis ===\n")
CKM_WORDS = ['aaB','AbA','AAb']

for idx in [1, 43]:
    M = snappy.OrientableClosedCensus[idx]
    coeffs, p, gens = get_h1_info(M)
    name = "m003 (PMNS)" if idx==1 else "m006 (CKM)"
    print(f"{name} idx={idx}: coeffs={list(coeffs)}, p={p}, gens={gens}")
    for w in CKM_WORDS:
        c = word_class(w, coeffs, gens, p)
        print(f"  [{w}] = {c}")
    ck = [word_class(w,coeffs,gens,p) for w in CKM_WORDS]
    print(f"  J_geom(aaB,AbA,AAb) = {j_geom(ck[0],ck[1],ck[2],p):.4f}")
    print(f"  Degenerate: {len(set(ck))<3}")
    print()

# Statistical summary
print("=== Census CP structure summary (indices 1-500) ===\n")
print(f"{'p':>5} {'count':>6} {'J_max_uniform':>14} {'MED_count':>10} {'J_min_med':>10}")
print("-"*50)

by_p = {}
import itertools

for idx in range(1, 501):
    try:
        M = snappy.OrientableClosedCensus[idx]
        coeffs, p, gens = get_h1_info(M)
        if coeffs is None: continue
        
        import itertools as it
        chars = list('abAB') if len(gens)==2 else list(gens[0]+gens[0].upper())
        words2 = [''.join(c) for c in it.product(chars[:4], repeat=2)][:20]
        words3 = [''.join(c) for c in it.product(chars[:4], repeat=3)][:30]
        all_w = [g for g in gens] + [g.upper() for g in gens] + words2 + words3
        
        classes = {}
        for w in all_w:
            try:
                classes[w] = word_class(w, coeffs, gens, p)
            except: pass
        
        wlist = list(classes.keys())
        j_max = 0.0
        for i in range(min(len(wlist),30)):
            for j in range(i+1,min(len(wlist),30)):
                for k in range(j+1,min(len(wlist),30)):
                    J = j_geom(classes[wlist[i]],classes[wlist[j]],
                               classes[wlist[k]],p)
                    if J > j_max: j_max = J
        
        if p not in by_p: by_p[p] = []
        by_p[p].append((idx, j_max))
    except: continue

# Theoretical max for small p
def theo_max(p):
    best = 0.0
    for k1 in range(min(p,10)):
        for k2 in range(k1+1,min(p,10)):
            for k3 in range(k2+1,min(p,10)):
                j = j_geom(k1,k2,k3,p)
                if j > best: best = j
    return best

for p in sorted(by_p.keys())[:25]:
    vals = by_p[p]
    j_t = theo_max(p)
    uniform = all(abs(v[1]-vals[0][1])<0.001 for v in vals)
    med = [v for v in vals if v[1] < j_t*0.95]
    print(f"{p:>5} {len(vals):>6} {'YES' if uniform else 'NO':>14} "
          f"{len(med):>10} {min(v[1] for v in med) if med else float('nan'):>10.3f}")
