import snappy
import numpy as np
import scipy.linalg as la
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore
import os

# ============================================================
# CONFIG
# ============================================================

OUT = "/mnt/c/dev/hyperbolic-flavor-scan/results"

phi = (1 + np.sqrt(5)) / 2
logphi = np.log(phi)

CUT = 6.0

SIGMAS = np.linspace(0.2, 1.5, 40)

MANIFOLDS = {
    "PMNS": ("m003", (-2,3)),
    "CKM":  ("m006", (-5,2))
}

# ============================================================
# UTILITIES
# ============================================================

def get_lengths(M, cutoff=CUT):
    spec = M.length_spectrum(cutoff)
    vals = []

    for s in spec:
        try:
            vals.append(float(s.length.real()))
        except:
            pass

    return np.array(sorted(vals))

def gaussian_kernel(L, sigma):
    D = squareform(pdist(L[:,None]))
    return np.exp(-(D**2)/(2*sigma*sigma))

def effective_rank(K):
    ev = np.linalg.eigvalsh(K)
    ev = np.maximum(ev, 1e-15)

    p = ev / ev.sum()

    return np.exp(-(p*np.log(p)).sum())

# ============================================================
# 1. LAPLACE-TYPE SPECTRUM
# ============================================================

def laplace_proxy(L, sigma=logphi):

    K = gaussian_kernel(L, sigma)

    D = np.diag(K.sum(axis=1))

    Lap = D - K

    ev = np.sort(np.linalg.eigvalsh(Lap))

    return ev

# ============================================================
# 2. DIRAC-LIKE OPERATOR
# ============================================================

def dirac_proxy(L):

    n = len(L)

    A = np.zeros((n,n))

    for i in range(n):
        for j in range(n):

            if i == j:
                continue

            d = abs(L[i]-L[j])

            A[i,j] = np.sin(d)/d if d > 1e-8 else 0

    A = (A + A.T)/2

    ev = np.sort(np.linalg.eigvalsh(A))

    return ev

# ============================================================
# 3. SELBERG-ZETA PROXY
# ============================================================

def selberg_proxy(L, s):

    total = 0.0 + 0.0j

    for ell in L:

        term = np.exp(-s*ell)

        total += np.log(abs(1-term))

    return np.real(total)

# ============================================================
# 4. TRANSFER-OPERATOR RESONANCES
# ============================================================

def transfer_operator(L, beta):

    n = len(L)

    T = np.zeros((n,n))

    for i in range(n):
        for j in range(n):

            d = abs(L[i]-L[j])

            T[i,j] = np.exp(-beta*d)

    return T

# ============================================================
# MAIN
# ============================================================

summary = []

for label, (name, fill) in MANIFOLDS.items():

    print("\n================================================")
    print(label)
    print("================================================")

    M = snappy.Manifold(name)
    M.dehn_fill(fill)

    L = get_lengths(M)

    print(f"Geodesics: {len(L)}")

    # --------------------------------------------------------
    # Laplace spectrum
    # --------------------------------------------------------

    lap_ev = laplace_proxy(L)

    print("\nLaplace proxy eigenvalues:")
    print(lap_ev[:10])

    # --------------------------------------------------------
    # Dirac spectrum
    # --------------------------------------------------------

    dir_ev = dirac_proxy(L)

    print("\nDirac proxy eigenvalues:")
    print(dir_ev[-10:])

    # --------------------------------------------------------
    # Effective-rank sweep
    # --------------------------------------------------------

    ranks = []

    for s in SIGMAS:

        K = gaussian_kernel(L, s)

        r = effective_rank(K)

        ranks.append(r)

    # --------------------------------------------------------
    # Selberg proxy
    # --------------------------------------------------------

    svals = np.linspace(0.2, 2.0, 200)

    zeta_vals = []

    for s in svals:
        zeta_vals.append(selberg_proxy(L, s))

    # --------------------------------------------------------
    # Transfer resonances
    # --------------------------------------------------------

    betas = np.linspace(0.2, 2.0, 60)

    radii = []

    for b in betas:

        T = transfer_operator(L, b)

        ev = np.linalg.eigvals(T)

        radii.append(np.max(np.abs(ev)))

    # --------------------------------------------------------
    # PLOTS
    # --------------------------------------------------------

    plt.figure(figsize=(8,5))
    plt.plot(SIGMAS, ranks)
    plt.axvline(logphi, linestyle='--')
    plt.xlabel("sigma")
    plt.ylabel("effective rank")
    plt.title(f"{label}: effective rank flow")
    plt.tight_layout()
    plt.savefig(f"{OUT}/{label}_rank_flow.png")
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(svals, zeta_vals)
    plt.xlabel("s")
    plt.ylabel("log|Z(s)| proxy")
    plt.title(f"{label}: Selberg proxy")
    plt.tight_layout()
    plt.savefig(f"{OUT}/{label}_selberg_proxy.png")
    plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(betas, radii)
    plt.xlabel("beta")
    plt.ylabel("spectral radius")
    plt.title(f"{label}: transfer resonances")
    plt.tight_layout()
    plt.savefig(f"{OUT}/{label}_transfer_resonance.png")
    plt.close()

    summary.append({
        "label": label,
        "rank_logphi":
            effective_rank(
                gaussian_kernel(L, logphi)
            ),
        "lap_gap":
            lap_ev[3]/lap_ev[2]
            if lap_ev[2] > 1e-10 else np.nan,
        "dirac_spread":
            np.std(dir_ev),
        "transfer_peak":
            max(radii)
    })

# ============================================================
# SUMMARY
# ============================================================

print("\n================================================")
print("SUMMARY")
print("================================================")

for s in summary:

    print(f"""
{s['label']}:

  effective rank @ log(phi):
      {s['rank_logphi']:.4f}

  Laplace spectral gap:
      {s['lap_gap']:.4f}

  Dirac spectral spread:
      {s['dirac_spread']:.4f}

  Max transfer resonance:
      {s['transfer_peak']:.4f}
""")

print("\nResults written to:")
print(OUT)

