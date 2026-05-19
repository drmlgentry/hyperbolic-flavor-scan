import snappy
import numpy as np

def get_coeffs_p(M):
    """Get homology coefficients for a,b generators and torsion order p."""
    h1 = M.homology()
    ed = h1.elementary_divisors()
    if len(ed) != 1:
        return None, None
    p = int(ed[0])
    if p < 3:
        return None, None
    
    fg = M.fundamental_group()
    gens = fg.generators()
    relators = fg.relators()
    n = len(gens)
    R = []
    for rel in relators:
        row = [0]*n
        for c in rel:
            idx = gens.index(c.lower())
            row[idx] += 1 if c.islower() else -1
        R.append(row)
    R = np.array(R, dtype=int)
    
    coeffs = [0]*n
    for i in range(n):
        if abs(R[0,i]) == p:
            coeffs[i] = 1
            for j in range(n):
                if j != i and R[1,j] != 0:
                    try:
                        inv = pow(int(R[1,j]), -1, p)
                        coeffs[j] = (-R[1,i] * inv) % p
                    except:
                        pass
            break
    return coeffs, p

def word_class(word, coeffs, p):
    val = 0
    for c in word:
        g = c.lower()
        if g == 'a': idx = 0
        elif g == 'b': idx = 1
        else: continue
        val += coeffs[idx] if c.islower() else -coeffs[idx]
    return val % p

def j_geom(k1, k2, k3, p):
    if k1==k2 or k2==k3 or k1==k3:
        return 0.0
    phi = [2*np.pi*k/p for k in [k1,k2,k3]]
    return abs(np.sin(phi[0]-phi[1]) +
               np.sin(phi[1]-phi[2]) +
               np.sin(phi[2]-phi[0]))

# Short words to use as geodesic candidates (length 1-3)
WORDS = ['a','b','A','B','aa','ab','aB','ba','bA','AB',
         'aab','aaB','abA','aBa','baa','baA','bAa']

print(f"{'Idx':>5} {'p':>4} {'Coeffs':>12} {'Triple':>20} "
      f"{'Classes':>15} {'J_geom':>8} {'Degenerate':>12}")
print("-"*85)

z5_results = []
for idx in [1, 2, 3, 43, 44, 45, 46, 47, 48]:  # m003, m006, and neighbors
    try:
        M = snappy.OrientableClosedCensus[idx]
        coeffs, p = get_coeffs_p(M)
        if coeffs is None:
            continue
        
        # Compute class for each word
        classes = {w: word_class(w, coeffs, p) for w in WORDS}
        
        # Find best 3-word combination by J_geom using actual geodesic classes
        best_J = -1
        best_combo = None
        words_list = list(classes.keys())
        for i in range(len(words_list)):
            for j in range(i+1, len(words_list)):
                for k in range(j+1, len(words_list)):
                    w1,w2,w3 = words_list[i],words_list[j],words_list[k]
                    k1,k2,k3 = classes[w1],classes[w2],classes[w3]
                    J = j_geom(k1,k2,k3,p)
                    if J > best_J:
                        best_J = J
                        best_combo = (w1,w2,w3,k1,k2,k3)
        
        deg = "YES" if best_combo and len({best_combo[3],best_combo[4],best_combo[5]})<3 else "no"
        print(f"{idx:>5} {p:>4} {str(coeffs):>12} "
              f"{str(best_combo[:3]) if best_combo else 'none':>20} "
              f"{str(best_combo[3:] if best_combo else ''):>15} "
              f"{best_J:>8.4f} {deg:>12}")
        
        if p == 5:
            z5_results.append((idx, p, coeffs, classes, best_J, best_combo))
            
    except Exception as e:
        print(f"{idx:>5}: error {e}")

# Show all Z/5 homology classes for m003 and m006
print("\n=== Z/5 word classes: m003 (idx 1) and m006 (idx 43) ===")
for idx, p, coeffs, classes, best_J, combo in z5_results:
    name = "m003" if idx==1 else f"idx{idx}"
    print(f"\n{name} (idx {idx}): coeffs={coeffs}, p={p}")
    print(f"  Word classes: {dict(list(classes.items())[:12])}")
    print(f"  Best J_geom triple: {combo}")
