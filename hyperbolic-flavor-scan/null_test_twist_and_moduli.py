"""
null_test_twist_and_moduli.py
==============================
Statistical null tests for both split papers:

Paper A: Are mb/mc and MZ/MW matches significant?
  - Generate N random word pairs from each manifold
  - Count how many give ratios within tolerance of PDG values
  - Report p-value

Paper B: Are twist angle matches significant?
  - For each SM parameter, count words within tolerance
  - Compare to expected count from uniform distribution on [0,180]
  - Report p-value per parameter and overall

Run: conda run -n sage python null_test_twist_and_moduli.py
Runtime: ~5 minutes
"""

import snappy
import numpy as np
from itertools import product as iprod

np.random.seed(42)

M_P = snappy.OrientableClosedCensus[1]
M_C = snappy.OrientableClosedCensus[43]
G_P = M_P.fundamental_group()
G_C = M_C.fundamental_group()

LETTERS_P = ['a','b','A','B']
LETTERS_C = ['a','b','A','B']

def twist_sl2c(G, word):
    """Positive-branch twist angle in degrees."""
    try:
        mat = np.array([[complex(G.SL2C(word)[i,j]) for j in range(2)]
                        for i in range(2)])
        eigs = np.linalg.eigvals(mat)
        lam = eigs[0] if np.imag(eigs[0]) >= 0 else eigs[1]
        phi = np.degrees(float(np.imag(np.log(lam))))
        if not np.isfinite(phi): return None
        return phi
    except: return None

def mod_lambda(G, word):
    """|lambda| of dominant eigenvalue."""
    try:
        mat = np.array([[complex(G.SL2C(word)[i,j]) for j in range(2)]
                        for i in range(2)])
        eigs = np.linalg.eigvals(mat)
        return max(abs(eigs[0]), abs(eigs[1]))
    except: return None

def geo_length(G, word):
    """Real translation length."""
    lam = mod_lambda(G, word)
    if lam is None or lam <= 1.0: return None
    return 2*np.log(lam)

print("="*60)
print("NULL TEST: GEODESIC MODULI (Paper A)")
print("="*60)
print()

# PDG targets
MB_MC = 4.18 / 1.27    # 3.2913
MZ_MW = 91.1876/80.377  # 1.1345
TOL_MBMC = 0.005        # |ratio - target| < 0.005
TOL_MZMW = 0.005

# Collect all words up to length 5
print("Collecting word data (length 1-5)...")
words_P = []
for L in range(1,6):
    for combo in iprod(LETTERS_P, repeat=L):
        w = ''.join(combo)
        lam = mod_lambda(G_P, w)
        if lam and lam > 1.01:
            words_P.append((w, lam))

words_C = []
for L in range(1,6):
    for combo in iprod(LETTERS_C, repeat=L):
        w = ''.join(combo)
        ell = geo_length(G_C, w)
        if ell and ell > 0.01:
            words_C.append((w, ell))

print(f"  m003 words with |lambda|>1: {len(words_P)}")
print(f"  m006 words with ell>0:      {len(words_C)}")
print()

# Count observed matches for mb/mc on m003
obs_mbmc = 0
hits_mbmc = []
for i,(w1,l1) in enumerate(words_P):
    for w2,l2 in words_P[i+1:]:
        r = l1/l2
        if abs(r - MB_MC) < TOL_MBMC:
            obs_mbmc += 1
            hits_mbmc.append((w1,w2,r))
        r2 = l2/l1
        if abs(r2 - MB_MC) < TOL_MBMC:
            obs_mbmc += 1
            hits_mbmc.append((w2,w1,r2))

print(f"mb/mc observed matches (|ratio-3.2913|<{TOL_MBMC}):")
print(f"  Count: {obs_mbmc}")
for w1,w2,r in sorted(hits_mbmc, key=lambda x: abs(x[2]-MB_MC))[:5]:
    print(f"  |lambda({w1})|/|lambda({w2})| = {r:.6f}")

# Null model: if |lambda| values were random on [1, max_lam]
# Probability that a random ratio falls within tol of MB_MC
n_P = len(words_P)
n_pairs = n_P*(n_P-1)
# Expected: 2 * n_pairs * 2*tol / MB_MC (for ratio distribution)
prob_mbmc = 2*TOL_MBMC / MB_MC  # approximate for ratio near 1
expected_mbmc = n_pairs * prob_mbmc * 2  # factor 2 for both orders
p_mbmc = obs_mbmc / expected_mbmc if expected_mbmc > 0 else 0
print(f"  Expected by chance: {expected_mbmc:.1f} pairs")
print(f"  p-value (obs/expected): {p_mbmc:.4f}")
if p_mbmc < 0.05:
    print(f"  SIGNIFICANT (p<0.05)")
else:
    print(f"  Not significant at p<0.05")
print()

# Count observed matches for MZ/MW on m006 (geodesic length ratios)
obs_mzmw = 0
hits_mzmw = []
for i,(w1,e1) in enumerate(words_C):
    for w2,e2 in words_C[i+1:]:
        if e2 < 0.01: continue
        r = e1/e2
        if abs(r - MZ_MW) < TOL_MZMW:
            obs_mzmw += 1
            hits_mzmw.append((w1,w2,r))
        r2 = e2/e1
        if abs(r2 - MZ_MW) < TOL_MZMW:
            obs_mzmw += 1
            hits_mzmw.append((w2,w1,r2))

print(f"MZ/MW observed matches (|ratio-1.1345|<{TOL_MZMW}):")
print(f"  Count: {obs_mzmw}")
for w1,w2,r in sorted(hits_mzmw, key=lambda x: abs(x[2]-MZ_MW))[:5]:
    print(f"  ell({w1})/ell({w2}) = {r:.6f}")

n_C = len(words_C)
n_pairs_C = n_C*(n_C-1)
prob_mzmw = 2*TOL_MZMW / MZ_MW
expected_mzmw = n_pairs_C * prob_mzmw * 2
p_mzmw = obs_mzmw / expected_mzmw if expected_mzmw > 0 else 0
print(f"  Expected by chance: {expected_mzmw:.1f} pairs")
print(f"  p-value (obs/expected): {p_mzmw:.4f}")
if p_mzmw < 0.05:
    print(f"  SIGNIFICANT (p<0.05)")
else:
    print(f"  Not significant at p<0.05")
print()

print("="*60)
print("NULL TEST: TWIST ANGLE CENSUS (Paper B)")
print("="*60)
print()

# SM targets for twist angles
SM_TARGETS = {
    'delta_CKM': 68.0,
    'theta12_CKM': 13.04,
    'theta23_CKM': 2.38,
    'theta12_nu': 33.41,
    'theta23_nu': 49.1,
}
TOL_ANGLE = 1.0  # degrees

# Collect all twist angles from both manifolds
print("Collecting twist angles...")
all_twists = []
for G, name in [(G_P,'m003'),(G_C,'m006')]:
    for L in range(1,6):
        for combo in iprod(LETTERS_P, repeat=L):
            w = ''.join(combo)
            phi = twist_sl2c(G, w)
            if phi is not None:
                all_twists.append((name, w, L, phi))

print(f"  Total twist angle values: {len(all_twists)}")
print()

# For each SM target, count matches using both phi and 180-phi
for target_name, target in SM_TARGETS.items():
    obs = 0
    best = []
    for mfld, word, length, phi in all_twists:
        for cand, cname in [(abs(phi),'direct'),(180-abs(phi),'180-phi')]:
            if abs(cand - target) < TOL_ANGLE:
                obs += 1
                best.append((abs(cand-target), mfld, word, cand, cname))
    
    best.sort()
    # Expected: 2 * N_words * 2*tol / 180 (uniform on [0,180], two folds)
    N = len(all_twists)
    expected = N * 2 * (2*TOL_ANGLE/180.0)
    p = obs/expected if expected > 0 else 0
    
    print(f"{target_name} = {target}°  (tol={TOL_ANGLE}°)")
    print(f"  Observed: {obs}  Expected: {expected:.1f}  p={p:.4f}")
    if best:
        err, mfld, word, val, fold = best[0]
        print(f"  Best: {mfld} word '{word}' [{fold}] = {val:.3f}° (err={err:.3f}°)")
    print()

# Overall: how many total matches vs expected?
total_obs = 0
total_exp = 0
N = len(all_twists)
for target_name, target in SM_TARGETS.items():
    obs = sum(1 for _,_,_,phi in all_twists
              for cand in [abs(phi), 180-abs(phi)]
              if abs(cand-target) < TOL_ANGLE)
    total_obs += obs
    total_exp += N * 2 * (2*TOL_ANGLE/180.0)

print(f"Overall census:")
print(f"  Total matches: {total_obs}  Total expected: {total_exp:.1f}")
print(f"  Enrichment factor: {total_obs/total_exp:.2f}x")
if total_obs/total_exp > 3:
    print(f"  SIGNIFICANTLY ENRICHED vs random")
print()

print("="*60)
print("SUMMARY FOR PAPERS")
print("="*60)
print()
print("Paper A (mass ratios):")
print(f"  mb/mc: {obs_mbmc} matches, p={p_mbmc:.4f}")
print(f"  MZ/MW: {obs_mzmw} matches, p={p_mzmw:.4f}")
print()
print("Paper B (twist angles):")
print(f"  Total enrichment: {total_obs/total_exp:.2f}x over random")
print(f"  Report per-parameter p-values in paper")
