"""
pmns_complex_scan.py
Complex Borel construction scan on m003 -- length-2 and length-3 words
Searches for word triples achieving simultaneous good fitness AND J near PDG

Target: fitness < 0.05, |J - J_pdg| < 0.005
J_pdg = -0.00966 (correct value at delta_CP = 197 deg)

Output: pmns_complex_scan_results.csv
"""

import snappy
import numpy as np
from scipy.linalg import logm, qr
from scipy.optimize import minimize
from itertools import permutations, product as iproduct
import pandas as pd
import time

# ── constants ────────────────────────────────────────────────────────
PMNS = np.array([
    [0.821, 0.550, 0.148],
    [0.357, 0.339, 0.871],
    [0.442, 0.762, 0.471]
])
J_PDG = -0.00966   # correct Jarlskog at delta_CP=197 deg
F_TARGET = 0.050   # fitness threshold to save result

# ── manifold setup ───────────────────────────────────────────────────
M = snappy.OrientableClosedCensus[1]   # m003
G = M.fundamental_group()
print("Manifold:", M.name(), " vol:", M.volume(), " H1:", M.homology())

# ── axis extraction ──────────────────────────────────────────────────
def get_complex_axis(word, alpha=1.0):
    """Extract complex axis from holonomy word.
    alpha=0: real only (Doc-8 method)
    alpha=1: full complex (Re + i*Im)
    """
    mat = np.array([[complex(G.SL2C(word)[i,j])
                     for j in range(2)] for i in range(2)])
    L = logm(mat)
    nr = np.array([
        np.real(L[0,1] + L[1,0]),
        np.real(1j*(L[1,0] - L[0,1])),
        np.real(L[0,0] - L[1,1])
    ])
    ni = np.array([
        np.imag(L[0,1] + L[1,0]),
        np.imag(1j*(L[1,0] - L[0,1])),
        np.imag(L[0,0] - L[1,1])
    ])
    nc = nr + alpha * 1j * ni
    nm = np.sqrt(np.real(np.dot(nc, np.conj(nc))))
    return nc / nm if nm > 1e-10 else nc

def evaluate(axes, scale):
    """Compute fitness and Jarlskog for given axes and scale."""
    l21 = scale * np.dot(np.conj(axes[0]), axes[1])
    l31 = scale * np.dot(np.conj(axes[0]), axes[2])
    l32 = scale * np.dot(np.conj(axes[1]), axes[2])
    L = np.array([[1+0j, 0, 0],
                  [l21,  1, 0],
                  [l31, l32, 1]], dtype=complex)
    try:
        Q, R = qr(L)
        d = np.diag(R)
        Q = Q * (d / np.abs(d))[np.newaxis, :]
        U_abs = np.abs(Q)
        best_f = np.inf
        best_p = None
        for p in permutations(range(3)):
            f = np.linalg.norm(U_abs[:, list(p)] - PMNS, 'fro')
            if f < best_f:
                best_f = f
                best_p = p
        Qp = Q[:, list(best_p)]
        J = np.imag(Qp[0,0]*Qp[1,1]*np.conj(Qp[0,1])*np.conj(Qp[1,0]))
        return best_f, J, Qp
    except Exception:
        return 99.0, 0.0, None

def optimize_scale(axes):
    """Find optimal scale minimizing joint loss."""
    def loss(p):
        f, J, _ = evaluate(axes, abs(p[0]))
        return f**2 + 15 * (J - J_PDG)**2
    best = (99.0, 0.0, 1.0)
    # coarse grid
    for sc in np.linspace(0.3, 6.0, 40):
        f, J, _ = evaluate(axes, sc)
        if f + 3*abs(J - J_PDG) < best[0] + 3*abs(best[1] - J_PDG):
            best = (f, J, sc)
    # Nelder-Mead polish
    try:
        res = minimize(loss, [best[2]], method='Nelder-Mead',
                       options={'maxiter': 1000, 'xatol': 1e-7, 'fatol': 1e-9})
        sc_opt = abs(res.x[0])
        f, J, Qp = evaluate(axes, sc_opt)
        return f, J, sc_opt, Qp
    except Exception:
        return best[0], best[1], best[2], None

# ── word lists ───────────────────────────────────────────────────────
words2 = [a+b for a in 'abAB' for b in 'abAB']
words3 = [a+b+c for a in 'abAB' for b in 'abAB' for c in 'abAB']

# Mixed triples: (len2, len2, len3) and (len2, len3, len3)
def gen_mixed_triples():
    # 2+2+3
    for w1, w2 in iproduct(words2, repeat=2):
        if w1 == w2: continue
        for w3 in words3:
            if w3 not in (w1, w2):
                yield (w1, w2, w3)
    # 2+3+3
    for w1 in words2:
        for w2, w3 in iproduct(words3, repeat=2):
            if len({w1, w2, w3}) == 3:
                yield (w1, w2, w3)

# Alpha values to try
ALPHAS = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

print("Starting complex Borel scan on m003...")
print(f"Targets: fitness < {F_TARGET}, |J - {J_PDG:.5f}| < 0.005")
print()

results = []
t0 = time.time()
n_scanned = 0
n_saved = 0

for triple in gen_mixed_triples():
    w1, w2, w3 = triple
    n_scanned += 1

    if n_scanned % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  Scanned {n_scanned} triples in {elapsed:.0f}s -- saved {n_saved}")

    for alpha in ALPHAS:
        try:
            axes = [get_complex_axis(w, alpha) for w in triple]
            f, J, sc, Qp = optimize_scale(axes)

            if f < F_TARGET:
                joint = f + 3*abs(J - J_PDG)
                results.append({
                    'w1': w1, 'w2': w2, 'w3': w3,
                    'alpha': round(alpha, 2),
                    'scale': round(sc, 5),
                    'fitness': round(f, 6),
                    'J': round(J, 6),
                    'J_err': round(abs(J - J_PDG), 6),
                    'joint_loss': round(joint, 6),
                })
                n_saved += 1

                # Print immediately if very good
                if f < 0.030 or abs(J - J_PDG) < 0.002:
                    print(f"*** GOOD: {w1}/{w2}/{w3} alpha={alpha} "
                          f"fit={f:.5f} J={J:.5f} scale={sc:.4f}")
        except Exception:
            pass

print()
print(f"Scan complete. {n_scanned} triples scanned, {n_saved} results saved.")

# Save results
df = pd.DataFrame(results)
if len(df) > 0:
    df = df.sort_values('joint_loss')
    outpath = r'C:\dev\hyperbolic-flavor-scan\pmns_complex_scan_results.csv'
    df.to_csv(outpath, index=False)
    print(f"Results saved to {outpath}")
    print()
    print("TOP 20 RESULTS:")
    print(df.head(20).to_string(index=False))
else:
    print("No results found within thresholds.")
    print("Consider raising F_TARGET or relaxing J tolerance.")
