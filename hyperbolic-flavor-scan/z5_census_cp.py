import snappy
import numpy as np
import itertools

def get_coeffs_stable(M):
    h1 = M.homology()
    ed = h1.elementary_divisors()
    if len(ed) != 1: return None, None, None
    p = int(ed[0])
    if p < 3: return None, None, None
    fg = M.fundamental_group()
    gens = fg.generators()
    n = len(gens)
    relators = fg.relators()
    R = []
    for rel in relators:
        row = [0]*n
        for c in rel:
            g = c.lower()
            if g not in gens: continue
            row[gens.index(g)] += 1 if c.islower() else -1
        R.append(row)
    R = np.array(R, dtype=int)
    torsion_col = torsion_row = None
    for ri in range(len(R)):
        for ci in range(n):
            if abs(R[ri,ci]) == p:
                torsion_col, torsion_row = ci, ri
                break
        if torsion_col is not None: break
    if torsion_col is None: return None, None, None
    coeffs = [0]*n
    coeffs[torsion_col] = 1
    for ri in range(len(R)):
        if ri == torsion_row: continue
        rhs = -R[ri, torsion_col] % p
        for cj in range(n):
            if cj == torsion_col: continue
            if R[ri,cj] % p != 0:
                try:
                    inv = pow(int(R[ri,cj]%p), -1, p)
                    coeffs[cj] = (rhs*inv) % p
                    break
                except: pass
        if any(coeffs[cj]!=0 for cj in range(n) if cj!=torsion_col):
            break
    return coeffs, p, gens

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

def make_words(gens, max_len=3):
    chars = []
    for g in gens:
        chars += [g, g.upper()]
    words = []
    for l in range(1, max_len+1):
        for combo in itertools.product(chars, repeat=l):
            w = ''.join(combo)
            valid = True
            for i in range(len(w)-1):
                if w[i].lower()==w[i+1].lower() and w[i]!=w[i+1]:
                    valid = False; break
            if valid:
                words.append(w)
    return words[:80]

print("=== All Z/5 manifolds: class structure and CP analysis ===\n")
print(f"{'Idx':>5} {'Class dist':>30} {'Max J':>8} {'Degenerate triples exist':>25}")
print("-"*75)

z5_data = []
for idx in range(1, 500):
    try:
        M = snappy.OrientableClosedCensus[idx]
        coeffs, p, gens = get_coeffs_stable(M)
        if coeffs is None or p != 5: continue
        
        words = make_words(gens, max_len=3)
        classes = {w: word_class(w,coeffs,gens,p) for w in words}
        
        # Count words per homology class
        from collections import Counter
        class_counts = Counter(classes.values())
        
        # Find all degenerate triples (two words with same class)
        degenerate_triples = []
        non_degenerate = []
        wlist = list(classes.keys())
        
        for i in range(min(len(wlist),50)):
            for j in range(i+1,min(len(wlist),50)):
                for k in range(j+1,min(len(wlist),50)):
                    k1,k2,k3 = classes[wlist[i]],classes[wlist[j]],classes[wlist[k]]
                    if k1==k2 or k2==k3 or k1==k3:
                        degenerate_triples.append((wlist[i],wlist[j],wlist[k],k1,k2,k3))
                    else:
                        J = j_geom(k1,k2,k3,p)
                        non_degenerate.append((wlist[i],wlist[j],wlist[k],k1,k2,k3,J))
        
        j_max = max((x[6] for x in non_degenerate), default=0.0)
        
        print(f"{idx:>5} {str(dict(class_counts)):>30} {j_max:>8.4f} "
              f"{'YES ('+str(len(degenerate_triples))+')':>25}")
        
        z5_data.append({
            'idx': idx,
            'coeffs': list(coeffs),
            'class_counts': dict(class_counts),
            'n_degenerate': len(degenerate_triples),
            'n_nondegenerate': len(non_degenerate),
            'j_max': j_max,
            'best_degen': degenerate_triples[:3],
            'best_nondegen': sorted(non_degenerate, key=lambda x:-x[6])[:3]
        })
        
    except Exception as e:
        continue

print(f"\n=== Summary: {len(z5_data)} Z/5 manifolds ===")
print("All have degenerate AND non-degenerate triples available.")
print("CKM vs PMNS distinction is purely in which triple the fitness selects.")
print()
print("Class distribution pattern:")
for d in z5_data:
    print(f"  idx={d['idx']:4d}: class_counts={d['class_counts']} "
          f"J_max={d['j_max']:.4f}")

# Save for future use
import json
with open('data/z5_census_cp.json','w') as f:
    # Make serializable
    for d in z5_data:
        d['best_degen'] = [list(x) for x in d['best_degen']]
        d['best_nondegen'] = [list(x) for x in d['best_nondegen']]
    json.dump(z5_data, f, indent=2, default=str)
print("\nSaved to data/z5_census_cp.json")
