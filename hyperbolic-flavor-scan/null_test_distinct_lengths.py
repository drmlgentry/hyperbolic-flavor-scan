"""
null_test_distinct_lengths.py
Proper null test for mb/mc and MZ/MW using DISTINCT geodesic lengths,
not words. Isospectral words are collapsed first.
Run: conda run -n sage python null_test_distinct_lengths.py
"""
import snappy
import numpy as np
from itertools import product as iprod

M_P = snappy.OrientableClosedCensus[1]   # m003
M_C = snappy.OrientableClosedCensus[43]  # m006
G_P = M_P.fundamental_group()
G_C = M_C.fundamental_group()
LETTERS = ['a','b','A','B']

def geo_length(G, word):
    try:
        mat = np.array([[complex(G.SL2C(word)[i,j]) for j in range(2)]
                        for i in range(2)])
        eigs = np.linalg.eigvals(mat)
        lam = max(abs(eigs[0]), abs(eigs[1]))
        if lam <= 1.001:
            return None
        return 2*np.log(lam)
    except:
        return None

def collect_distinct(G, max_len=5):
    """Return list of (length, shortest_word) for distinct geodesic lengths."""
    seen = {}
    for L in range(1, max_len+1):
        for combo in iprod(LETTERS, repeat=L):
            w = ''.join(combo)
            ell = geo_length(G, w)
            if ell is None or ell < 0.01:
                continue
            key = round(ell, 4)
            if key not in seen:
                seen[key] = (ell, w)
    return sorted(seen.values())

def null_test_ratio(distinct_lengths, target, tol, label):
    """Test if any pair of distinct lengths gives ratio near target."""
    hits = []
    n = len(distinct_lengths)
    for i in range(n):
        e1, w1 = distinct_lengths[i]
        for j in range(i+1, n):
            e2, w2 = distinct_lengths[j]
            r = e1/e2 if e1 > e2 else e2/e1
            wnum = w1 if e1 > e2 else w2
            wden = w2 if e1 > e2 else w1
            if abs(r - target) < tol:
                hits.append((abs(r-target), wnum, wden, r))

    n_pairs = n*(n-1)//2
    # Null model: ratios uniform near 1. For ratio r = e1/e2 ~ target,
    # the probability a random pair lands within tol is approximately:
    # P = 2*tol / target  (for ratios distributed near target)
    expected = n_pairs * 2*tol / target
    p = len(hits)/expected if expected > 0 else 0

    hits.sort()
    print(f"\n{label} = {target:.5f}  (tol={tol})")
    print(f"  Distinct length pairs: {n_pairs}")
    print(f"  Observed hits: {len(hits)}")
    print(f"  Expected by chance: {expected:.2f}")
    print(f"  Enrichment: {len(hits)/expected:.2f}x")
    if expected > 0:
        # Binomial p-value (one-sided)
        from scipy.stats import binom
        p_binom = 1 - binom.cdf(len(hits)-1, n_pairs, 2*tol/target)
        print(f"  Binomial p-value: {p_binom:.4f}")
    print(f"  Top matches:")
    for err, wn, wd, r in hits[:5]:
        print(f"    ell({wn})/ell({wd}) = {r:.6f}  err={err:.6f}")
    return len(hits), expected


print("="*60)
print("DISTINCT GEODESIC LENGTH NULL TESTS")
print("="*60)

print("\nCollecting distinct lengths on m003...")
dist_P = collect_distinct(G_P, max_len=5)
print(f"  {len(dist_P)} distinct lengths from words of length 1-5")

print("\nCollecting distinct lengths on m006...")
dist_C = collect_distinct(G_C, max_len=5)
print(f"  {len(dist_C)} distinct lengths from words of length 1-5")

print("\n" + "="*60)
print("TEST 1: mb/mc on m003")
print("="*60)
MB_MC = 4.18/1.27
null_test_ratio(dist_P, MB_MC, tol=0.005, label="mb/mc")

print("\n" + "="*60)
print("TEST 2: MZ/MW on m006")
print("="*60)
MZ_MW = 91.1876/80.377
null_test_ratio(dist_C, MZ_MW, tol=0.005, label="MZ/MW")

print("\n" + "="*60)
print("TEST 3: mb/mc on m006 (cross-check)")
print("="*60)
null_test_ratio(dist_C, MB_MC, tol=0.005, label="mb/mc cross")

print("\n" + "="*60)
print("TEST 4: Tighter tolerance for best hits")
print("="*60)
# Use tol=0.001 (3x tighter) to see if best hits survive
null_test_ratio(dist_P, MB_MC, tol=0.001, label="mb/mc tight")
null_test_ratio(dist_C, MZ_MW, tol=0.001, label="MZ/MW tight")

print("\n" + "="*60)
print("SPECTRUM LISTING (first 20 distinct lengths, m003)")
print("="*60)
print(f"\n{'ell':>10}  {'word':>12}")
for ell, w in dist_P[:20]:
    print(f"  {ell:.6f}    {w}")

print("\n" + "="*60)
print("SPECTRUM LISTING (first 20 distinct lengths, m006)")
print("="*60)
print(f"\n{'ell':>10}  {'word':>12}")
for ell, w in dist_C[:20]:
    print(f"  {ell:.6f}    {w}")
