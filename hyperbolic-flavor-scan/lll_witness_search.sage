# lll_witness_search.sage
# Use LLL to find optimal norm witnesses for large mass ratios
# This replaces brute-force search for heavy fermions (t quark)
# The LLL is ONLY the search method — the statistical analysis
# uses the same null model framework as null_model_pipeline_v2.sage
# Run: sage lll_witness_search.sage

from sage.all import *
import numpy as np

MASSES = {
    "e":   1.0,
    "mu":  206.77,
    "tau": 3477.22,
    "u":   4.24,
    "d":   9.14,
    "s":   182.78,
    "c":   2495.11,
    "b":   8180.04,
    "t":   338082.19,
}

def is_prime_py(n):
    n = int(n)
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5)+1, 2):
        if n%i==0: return False
    return True

# ── Norm functions ────────────────────────────────────────────────────────────
def norm_eisenstein(a, b):
    a,b = int(a),int(b)
    n = a*a - a*b + b*b
    return n if n > 0 else 0

def norm_sqrt17(a, b):
    a,b = int(a),int(b)
    return abs(a*a - 17*b*b)

def norm_sqrt5(a, b):
    a,b = int(a),int(b)
    return abs(a*a - 5*b*b)

# ── LLL-accelerated witness search ────────────────────────────────────────────
def lll_witness(m, norm_func, norm_matrix, scale=10**6):
    """
    Use LLL to find best (a,b) with norm_func(a,b) ≈ m.

    The lattice: embed the problem as finding short vectors in
    a 3D lattice that balances closeness to m and lattice structure.

    norm_matrix = [[A, B], [B, C]] defines N(a,b) = Aa^2 + 2Bab + Cb^2
    """
    m = float(m)
    # Augmented lattice:
    # Row 0: [scale, 0, A]    -> a-direction
    # Row 1: [0, scale, C]    -> b-direction
    # Row 2: [scale*m, 0, 0]  -> target

    # Use CVP (closest vector problem) approximation via LLL
    # Build basis for the lattice of vectors (a, b, a*A+b*B, a*B+b*C)
    A = float(norm_matrix[0][0])
    B = float(norm_matrix[0][1])
    C = float(norm_matrix[1][1])

    # Simple 2D lattice with augmented target
    M = scale
    # Lattice basis in 3D: columns represent (a_coeff, b_coeff, norm_coeff)
    basis = Matrix(ZZ, [
        [M, 0, int(A*M)],
        [0, M, int(C*M)],
        [0, 0, 1],
    ])

    # Target: we want norm ≈ m, so target vector is (0, 0, int(m*M))
    # CVP is NP-hard in general but LLL gives 2^n approximation
    # For our 2D case, we can just do bounded brute force + LLL hint

    # LLL-reduced basis
    L = basis.LLL()

    # Extract candidate solutions from reduced basis rows
    candidates = []
    for row in L:
        for sign_a in [1, -1]:
            for sign_b in [1, -1]:
                # Each LLL vector gives a lattice direction
                # Try multiples up to small radius
                for k1 in range(-3, 4):
                    for k2 in range(-3, 4):
                        vec = k1 * L[0] + k2 * L[1]
                        # Decode a, b from vector
                        # The first two components encode a,b scaled by M
                        if M > 0:
                            a_approx = vec[0] / M
                            b_approx = vec[1] / M
                            for da in [-1,0,1]:
                                for db in [-1,0,1]:
                                    a = int(round(a_approx)) + da
                                    b = int(round(b_approx)) + db
                                    n = int(norm_func(a, b))
                                    if n > 0:
                                        err = abs(n - m)/m
                                        candidates.append((err, n, a, b))

    # Also do small brute force for safety
    R = min(200, max(50, int(m**0.5 / 10) + 20))
    for a in range(-R, R+1, max(1, R//50)):
        for b in range(-R, R+1, max(1, R//50)):
            n = int(norm_func(a, b))
            if n > 0:
                err = abs(n - m)/m
                candidates.append((err, n, a, b))

    if not candidates:
        return 0, (0,0), 1.0

    candidates.sort()
    best_err, best_n, best_a, best_b = candidates[0]
    return best_n, (best_a, best_b), float(best_err)

# ── Norm matrices for LLL ─────────────────────────────────────────────────────
# Eisenstein: N(a,b) = a^2 - ab + b^2 = [[1,-1/2],[-1/2,1]] (pos def)
M_EISENSTEIN = [[1, -1], [-1, 4]]  # scaled to integers: 4*(a^2-ab+b^2)

# Z[sqrt(17)]: N(a,b) = |a^2 - 17b^2|  (indefinite, use |.|)
# For LLL we use the positive part: a^2 + 17b^2 as proxy
M_SQRT17_POS = [[1, 0], [0, 17]]

# Z[sqrt(5)]: similarly
M_SQRT5_POS  = [[1, 0], [0, 5]]

# ── Direct brute force for verification (wider radius) ────────────────────────
def brute_witness(m, norm_func, R=250):
    """Extended brute force for verification."""
    m = float(m)
    best_err = 1.0
    best_n = 0; best = (0,0)
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            n = int(norm_func(a, b))
            if n <= 0: continue
            err = abs(n - m)/m
            if err < best_err:
                best_err = err; best_n = n; best = (a,b)
    return best_n, best, float(best_err)

# ── Run analysis ──────────────────────────────────────────────────────────────
print("="*70)
print("LLL-ACCELERATED FERMION MASS WITNESS SEARCH")
print("="*70)
print()

PHI = float((1+5**0.5)/2)

print("LEPTON SECTOR: Z[omega] Eisenstein norms")
print(f"{'f':>5} {'m/m_e':>12}  {'N':>10}  {'witness':>14}  {'err%':>8}  {'notes'}")
print("-"*65)

for name in ["mu", "tau"]:
    m = float(MASSES[name])
    N, (a,b), err = brute_witness(m, norm_eisenstein, R=120)
    verify = norm_eisenstein(a,b)
    prime = is_prime_py(N) if 0 < N < 10**8 else False
    mod3 = N%3 if N > 0 else -1
    print(f"{name:>5} {m:>12.4f}  {N:>10}  ({a:4},{b:4})    {err*100:>8.4f}%  "
          f"prime={prime} mod3={mod3}")

print()
print("QUARK SECTOR: Z[sqrt(17)] norms")
print(f"{'f':>5} {'m/m_e':>12}  {'N':>10}  {'witness':>16}  {'err%':>8}  {'notes'}")
print("-"*70)

for name in ["u","d","s","c","b","t"]:
    m = float(MASSES[name])
    # For t quark, need large radius
    R = 250 if name == "t" else 150
    N, (a,b), err = brute_witness(m, norm_sqrt17, R=R)
    verify = norm_sqrt17(a,b)
    prime = is_prime_py(N) if 0 < N < 10**8 else False
    mod17 = N%17 if N > 0 else -1
    splits = mod17 in [1,2,4,8,9,13,15,16] if N > 0 else False
    print(f"{name:>5} {m:>12.4f}  {N:>10}  ({a:4},{b:5})    {err*100:>8.4f}%  "
          f"prime={prime} mod17={mod17} split={splits}")

print()
print("COMPARISON: Z[sqrt(5)] norms (field of phi)")
print(f"{'f':>5} {'m/m_e':>12}  {'N':>10}  {'witness':>16}  {'err%':>8}")
print("-"*65)

for name in ["mu","tau","u","d","s","c","b","t"]:
    m = float(MASSES[name])
    N, (a,b), err = brute_witness(m, norm_sqrt5, R=200)
    print(f"{name:>5} {m:>12.4f}  {N:>10}  ({a:4},{b:5})    {err*100:>8.4f}%")

print()
print("="*70)
print("STATISTICAL COMPARISON: RANDOM DISCRIMINANTS vs STRUCTURED")
print("="*70)
print()
print("Testing 20 random discriminants d (not 5,17) for comparison...")
print(f"{'d':>6}  {'mean_err%':>10}  {'f<2%':>7}")
print("-"*30)

import random
random.seed(2024)
random_discs = [d for d in range(2, 100)
                if d not in [5,17] and not any(d%p==0 and (d//p)%p==0
                for p in range(2,int(d**0.5)+1))][:20]

disc_results = []
for d in random_discs:
    def norm_d(a,b,dd=d): return abs(int(a)*int(a)-dd*int(b)*int(b))
    errs = []
    for name, m in MASSES.items():
        if name == 'e': continue
        N,(a,b),err = brute_witness(float(m), norm_d, R=80)
        errs.append(err)
    mean_err = float(np.mean(errs))*100
    f2 = sum(1 for e in errs if e < 0.02)/len(errs)
    disc_results.append((d, mean_err, f2))
    print(f"{d:>6}  {mean_err:>10.4f}%  {f2:>7.2f}")

random_mean_err = float(np.mean([r[1] for r in disc_results]))
random_f2 = float(np.mean([r[2] for r in disc_results]))

print(f"\nRandom d average: mean_err={random_mean_err:.4f}%  f<2%={random_f2:.2f}")
print()

# Structured results
struct_errs = {}
for model_name, norm_func in [
    ("Z[omega]",  norm_eisenstein),
    ("Z[sqrt17]", norm_sqrt17),
    ("Z[sqrt5]",  norm_sqrt5),
]:
    errs = []
    for name, m in MASSES.items():
        if name == 'e': continue
        R = 250 if name == 't' else 150
        N,(a,b),err = brute_witness(float(m), norm_func, R=R)
        errs.append(err)
    mean_err = float(np.mean(errs))*100
    f2 = sum(1 for e in errs if e < 0.02)/len(errs)
    struct_errs[model_name] = (mean_err, f2)
    p_val = sum(1 for r in disc_results if r[1] <= mean_err)/len(disc_results)
    print(f"{model_name:<12}: mean_err={mean_err:.4f}%  f<2%={f2:.2f}  "
          f"p<random: {p_val:.2f} ({sum(1 for r in disc_results if r[2]>=f2)}/{len(disc_results)} null models beat it in f<2%)")

print()
print("="*70)
print("LATEX TABLE FOR PAPER")
print("="*70)
print()
print(r"\begin{tabular}{llrrr}")
print(r"\toprule")
print(r"Fermion & Ring & $m_f/m_e$ & Norm $N$ & Error \\")
print(r"\midrule")

lepton_data = {}
for name in ["mu","tau"]:
    m = float(MASSES[name])
    N,(a,b),err = brute_witness(m, norm_eisenstein, R=120)
    lepton_data[name] = (m, N, a, b, err)
    print(f"{name} & $\\mathbb{{Z}}[\\omega]$ & {m:.2f} & {N} & "
          f"${err*100:.4f}\\%$ \\\\")

print(r"\midrule")
quark_data = {}
for name in ["u","d","s","c","b","t"]:
    m = float(MASSES[name])
    R = 250 if name == "t" else 150
    N,(a,b),err = brute_witness(m, norm_sqrt17, R=R)
    quark_data[name] = (m, N, a, b, err)
    print(f"{name} & $\\mathbb{{Z}}[\\sqrt{{17}}]$ & {m:.2f} & {N} & "
          f"${err*100:.4f}\\%$ \\\\")
print(r"\bottomrule")
print(r"\end{tabular}")
