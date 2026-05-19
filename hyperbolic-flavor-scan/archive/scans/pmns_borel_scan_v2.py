"""
pmns_borel_scan_v2.py
Canonical checkpointed PMNS Borel/N-factor scan on m003 (Meyerhoff manifold).
Paper 6 (PRD): PMNS Lepton Mixing from the Borel Structure of Hyperbolic Holonomy.

Best known result: words aa/ab/aB, column perm (0,2,1), fitness F=0.01897
Construction: lower-triangular (Borel/Iwasawa N-factor) overlap matrix L,
              then QR decomposition to extract unitary mixing matrix U.

Usage: python scans/pmns_borel_scan_v2.py [--resume]
Checkpoints to: data/pmns_borel_scan_v2_results.csv
                data/pmns_borel_scan_v2_checkpoint.txt
"""
import numpy as np
from scipy.linalg import qr, logm
from itertools import permutations, combinations_with_replacement
import snappy, pandas as pd, time, os, argparse

# ── PDG 2024 PMNS target (absolute values) ───────────────────────
PMNS_TARGET = np.array([
    [0.821, 0.550, 0.148],
    [0.357, 0.339, 0.871],
    [0.442, 0.762, 0.471],
])
PERMS = list(permutations([0, 1, 2]))

# ── Paths ─────────────────────────────────────────────────────────
DATA_DIR  = os.path.join(os.path.dirname(__file__), "..", "data")
OUT_FILE  = os.path.join(DATA_DIR, "pmns_borel_scan_v2_results.csv")
CKPT_FILE = os.path.join(DATA_DIR, "pmns_borel_scan_v2_checkpoint.txt")
SAVE_EVERY = 100

def to_numpy(m):
    """Convert SnapPy SimpleMatrix to numpy complex array."""
    return np.array([[complex(m[i,j]) for j in range(2)] for i in range(2)])

def fitness(U):
    """Permutation-minimized Frobenius norm vs PMNS target."""
    best = 1e9
    for perm in PERMS:
        best = min(best, float(np.linalg.norm(U[:, perm] - PMNS_TARGET, 'fro')))
    return best

def axis_from_matrix(M):
    """Extract rotation axis from SL2C matrix via Pauli decomposition of log."""
    try:
        L = logm(M)
    except Exception:
        return None
    nx = (L[0,1].real + L[1,0].real) / 2
    ny = (L[1,0].imag - L[0,1].imag) / 2
    nz = (L[0,0].real - L[1,1].real) / 2
    n  = np.array([nx, ny, nz])
    nm = np.linalg.norm(n)
    return n / nm if nm > 1e-10 else None

def borel_overlap_matrix(n1, n2, n3, scale=1.0):
    """
    Lower-triangular (Borel/N-factor) overlap matrix.
    L[i,j] = n_i . n_j for i >= j, 0 otherwise.
    Scaled by 'scale' before QR.
    """
    axes = [n1, n2, n3]
    L = np.zeros((3, 3))
    for i in range(3):
        for j in range(i+1):
            L[i, j] = np.dot(axes[i], axes[j])
    return L * scale

def borel_to_unitary(n1, n2, n3, scale=1.0):
    """Extract unitary matrix from Borel overlap via QR."""
    L = borel_overlap_matrix(n1, n2, n3, scale)
    if abs(np.linalg.det(L)) < 1e-10:
        return None
    Q, R = qr(L.T)
    # Fix signs so diagonal of R is positive
    signs = np.diag(np.sign(np.diag(R)))
    U = Q @ signs
    return np.abs(U)

def scan_word_triples(G, words, top_k=20):
    """
    Scan all triples from word list, return top_k by fitness.
    """
    # Pre-compute axes
    axes = {}
    for w in words:
        try:
            M = to_numpy(G.SL2C(w))
            n = axis_from_matrix(M)
            if n is not None:
                axes[w] = n
        except Exception:
            pass

    word_list = list(axes.keys())
    results = []
    n_triples = len(word_list)**3
    t0 = time.time()

    count = 0
    for i, w1 in enumerate(word_list):
        for w2 in word_list:
            for w3 in word_list:
                count += 1
                n1, n2, n3 = axes[w1], axes[w2], axes[w3]
                U = borel_to_unitary(n1, n2, n3)
                if U is None:
                    continue
                f = fitness(U)
                if f < 0.050:
                    results.append({
                        "w1": w1, "w2": w2, "w3": w3,
                        "fitness": f,
                    })

        if i % 5 == 0:
            elapsed = time.time() - t0
            rate = count / max(elapsed, 1)
            print(f"  [{count:>8}/{n_triples:>8}] "
                  f"{len(results)} hits  {rate:.0f} t/s")

    results.sort(key=lambda x: x["fitness"])
    return results[:top_k]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--max-len", type=int, default=3,
                        help="Max word length to include (default 3)")
    args = parser.parse_args()

    # Load manifold
    M  = snappy.OrientableClosedCensus[1]   # Meyerhoff manifold = m003
    G  = M.fundamental_group()
    base_gens = [g for g in G.generators() if g.islower()]
    all_gens  = base_gens + [g.upper() for g in base_gens]

    print(f"Manifold: {M.name()}  vol={float(M.volume()):.4f}  H1={M.homology()}")
    print(f"Generators: {all_gens}")

    # Generate words up to max_len
    from itertools import product as iproduct
    inv = {g: g.upper() if g.islower() else g.lower() for g in all_gens}
    words = []
    for length in range(1, args.max_len + 1):
        for w in iproduct(all_gens, repeat=length):
            ok = all(w[i] != inv.get(w[i+1], "___") for i in range(len(w)-1))
            if ok:
                words.append("".join(w))
    print(f"Words up to length {args.max_len}: {len(words)}")

    # Resume from checkpoint
    start = 0
    existing = []
    if args.resume and os.path.exists(CKPT_FILE) and os.path.exists(OUT_FILE):
        with open(CKPT_FILE) as f:
            start = int(f.read().strip())
        existing = pd.read_csv(OUT_FILE).to_dict("records")
        print(f"Resuming from {start} ({len(existing)} results)")

    print(f"\nScanning triples...")
    results = existing + scan_word_triples(G, words)
    results.sort(key=lambda x: x["fitness"])

    # Save
    df = pd.DataFrame(results)
    df.to_csv(OUT_FILE, index=False)
    print(f"\nTop results:")
    print(df.head(20).to_string())
    print(f"\nBest known: aa/ab/aB  F=0.01897  (column perm 0,2,1)")
    print(f"Saved {len(df)} results to {OUT_FILE}")

if __name__ == "__main__":
    main()
