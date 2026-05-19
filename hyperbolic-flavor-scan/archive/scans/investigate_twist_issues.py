"""
investigate_twist_issues.py
Resolve three open issues from verify_twist_paper.py:
1. Old CP formula -- find correct sign combination for 203.5
2. MZ/MW ratio -- find correct formula
3. PMNS theta13 -- find correct folding convention
Run: conda run -n sage python investigate_twist_issues.py
"""
import snappy, numpy as np

M_P = snappy.OrientableClosedCensus[1]
M_C = snappy.OrientableClosedCensus[43]
rho_P = M_P.polished_holonomy()
rho_C = M_C.polished_holonomy()

def eigendata(rho, word):
    mat = np.array(rho(word), dtype=complex)
    mat /= np.sqrt(np.linalg.det(mat))
    evals = np.linalg.eigvals(mat)
    lam = evals[np.argmax(np.abs(evals))]
    ell = 2*abs(float(np.real(np.log(lam))))
    phi = np.degrees(float(np.imag(np.log(lam))))
    return abs(lam), ell, phi

# ── ISSUE 1: Old CP formula ───────────────────────────────────────
print("="*60)
print("ISSUE 1: Old CP formula -- find sign combo giving 203.5")
print("="*60)
words_cp = ['aa','ab','aB']
phis_cp = []
for w in words_cp:
    _, _, phi = eigendata(rho_P, w)
    phis_cp.append(phi)
    print(f"  phi({w} on m003) = {phi:.4f}")

print()
print("All sign combinations (mod 360):")
for s1 in [1,-1]:
    for s2 in [1,-1]:
        for s3 in [1,-1]:
            for offset in [0, 180, -180, 360, -360]:
                val = (s1*phis_cp[0] + s2*phis_cp[1] + s3*phis_cp[2] + offset) % 360
                if abs(val-203.5) < 2.0:
                    print(f"  {s1}*phi(aa)+{s2}*phi(ab)+{s3}*phi(aB)+{offset} "
                          f"= {val:.4f} <-- MATCH 203.5")

# Also try the other eigenvalue (conjugate)
print()
print("Trying conjugate eigenvalues (phi -> -phi):")
for s1 in [1,-1]:
    for s2 in [1,-1]:
        for s3 in [1,-1]:
            for offset in [0, 180, 360]:
                val = (s1*(-phis_cp[0]) + s2*(-phis_cp[1]) + s3*(-phis_cp[2]) + offset) % 360
                if abs(val-203.5) < 2.0:
                    print(f"  {s1}*(-phi(aa))+{s2}*(-phi(ab))+{s3}*(-phi(aB))+{offset} "
                          f"= {val:.4f} <-- MATCH 203.5")

# ── ISSUE 2: MZ/MW ratio ─────────────────────────────────────────
print()
print("="*60)
print("ISSUE 2: MZ/MW ratio -- find correct formula")
print("="*60)
MZ_MW = 91.1876 / 80.377  # PDG 2024 = 1.13451

print(f"\nPDG MZ/MW = {MZ_MW:.6f}")
print()

# Check eigenvalue magnitudes and lengths for many word pairs
print("Scanning word pairs for |lambda| ratio near MZ/MW:")
from itertools import product as iprod

letters = ['a','b','A','B']
words2 = [''.join(p) for p in iprod(letters, repeat=2)]
words3 = [''.join(p) for p in iprod(letters, repeat=3)]
words4 = [''.join(p) for p in iprod(letters, repeat=4)]
words5 = [''.join(p) for p in iprod(letters, repeat=5)]
all_words = words2 + words3 + words4 + words5

# Get data for all words on both manifolds
data_C = {}
for w in all_words:
    try:
        lam, ell, phi = eigendata(rho_C, w)
        if ell > 0.01:  # skip identity-like
            data_C[w] = (lam, ell, phi)
    except:
        pass

print(f"  Got data for {len(data_C)} words on m006")
print()

# Check ratio of |lambda| for all pairs
hits = []
for w1, (l1, e1, p1) in data_C.items():
    for w2, (l2, e2, p2) in data_C.items():
        if w1 >= w2: continue
        if l2 < 0.01: continue
        r = l1/l2
        if abs(r - MZ_MW) < 0.01:
            hits.append((abs(r-MZ_MW), w1, w2, r, 'lam_ratio'))
        r2 = np.exp(e1)/np.exp(e2)
        if abs(r2 - MZ_MW) < 0.01:
            hits.append((abs(r2-MZ_MW), w1, w2, r2, 'exp_ell_ratio'))
        r3 = e1/e2
        if abs(r3 - MZ_MW) < 0.01:
            hits.append((abs(r3-MZ_MW), w1, w2, r3, 'ell_ratio'))

hits.sort()
if hits:
    print("Top matches for MZ/MW on m006:")
    for err, w1, w2, r, method in hits[:10]:
        print(f"  {method}: {w1}/{w2} = {r:.6f} (err={err:.6f})")
else:
    print("No matches found for MZ/MW on m006 with word length ≤ 5")

# Also check m003
data_P = {}
for w in all_words:
    try:
        lam, ell, phi = eigendata(rho_P, w)
        if ell > 0.01:
            data_P[w] = (lam, ell, phi)
    except:
        pass

hits_P = []
for w1, (l1, e1, p1) in data_P.items():
    for w2, (l2, e2, p2) in data_P.items():
        if w1 >= w2: continue
        if l2 < 0.01: continue
        for r, method in [(l1/l2,'lam'),(e1/e2,'ell')]:
            if abs(r - MZ_MW) < 0.01:
                hits_P.append((abs(r-MZ_MW), w1, w2, r, method))
hits_P.sort()
if hits_P:
    print("\nTop matches for MZ/MW on m003:")
    for err, w1, w2, r, method in hits_P[:5]:
        print(f"  {method}: {w1}/{w2} = {r:.6f} (err={err:.6f})")
else:
    print("No matches on m003 either")

# ── ISSUE 3: PMNS theta13 ─────────────────────────────────────────
print()
print("="*60)
print("ISSUE 3: PMNS theta13=8.57 -- find correct formula")
print("="*60)
PDG_t13 = 8.57  # degrees

print("\nTrying all words up to length 3 on m003:")
words3_P = [''.join(p) for p in iprod(letters, repeat=1)]
words3_P += [''.join(p) for p in iprod(letters, repeat=2)]
words3_P += [''.join(p) for p in iprod(letters, repeat=3)]

hits_t13 = []
for w in words3_P:
    try:
        lam, ell, phi = eigendata(rho_P, w)
        if ell < 0.01: continue
        # Try all folding conventions
        candidates = {
            '|phi|': abs(phi),
            '180-|phi|': 180-abs(phi),
            'phi%180': phi%180,
            '|phi|%90': abs(phi)%90,
            '(180+phi)%360': (180+phi)%360,
            'ell*r2d': ell*180/np.pi,
        }
        for name, val in candidates.items():
            if abs(val - PDG_t13) < 0.5:
                hits_t13.append((abs(val-PDG_t13), w, name, val))
    except:
        pass

hits_t13.sort()
print(f"\nBest matches for theta13={PDG_t13}° on m003:")
for err, w, formula, val in hits_t13[:15]:
    print(f"  {w}: {formula} = {val:.4f}° (err={err:.4f}°)")

# Also check m006
print("\nBest matches for theta13 on m006:")
hits_t13_C = []
words_scan = [''.join(p) for p in iprod(letters, repeat=1)]
words_scan += [''.join(p) for p in iprod(letters, repeat=2)]
words_scan += [''.join(p) for p in iprod(letters, repeat=3)]
for w in words_scan:
    try:
        lam, ell, phi = eigendata(rho_C, w)
        if ell < 0.01: continue
        for name, val in [('|phi|',abs(phi)),('180-|phi|',180-abs(phi))]:
            if abs(val - PDG_t13) < 0.5:
                hits_t13_C.append((abs(val-PDG_t13), w, name, val))
    except:
        pass
hits_t13_C.sort()
for err, w, formula, val in hits_t13_C[:10]:
    print(f"  {w}: {formula} = {val:.4f}° (err={err:.4f}°)")

print()
print("="*60)
print("SUMMARY")
print("="*60)
print("Results above determine:")
print("1. Whether 203.5 is reproducible (and with what formula)")
print("2. Whether MZ/MW has any geometric origin in these manifolds")
print("3. Which word and formula gives theta13=8.57")
