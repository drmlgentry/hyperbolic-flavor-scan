import snappy
import numpy as np
from scipy.linalg import logm, qr
from itertools import product, combinations
import pandas as pd
import time
import os

# ── PMNS target (PDG 2024) ────────────────────────────────────────
PMNS = np.array([
    [0.821, 0.550, 0.148],
    [0.357, 0.339, 0.871],
    [0.442, 0.762, 0.471]
])

CENSUS_RANGE  = 200
MAX_WORDLEN   = 4
SIGMA_VALUES  = np.linspace(0.30, 1.50, 25)
TOP_CANDIDATES = 500
TOP_N         = 50
RESULTS_FILE  = "pmns_scan_results.csv"
PROGRESS_FILE = "pmns_scan_progress.txt"

# ── Word generation ───────────────────────────────────────────────
def gen_words(max_len):
    letters = ['a','A','b','B']
    cancel  = {'a':'A','A':'a','b':'B','B':'b'}
    words   = []
    for length in range(2, max_len+1):
        for w in product(letters, repeat=length):
            bad = any(cancel.get(w[i]) == w[i+1] for i in range(len(w)-1))
            if not bad:
                words.append(''.join(w))
    return words

# ── Axis extraction ───────────────────────────────────────────────
def matrix_to_axis(mat):
    mat = np.array(mat, dtype=complex)
    det = np.linalg.det(mat)
    if abs(det) < 1e-12:
        return None
    mat = mat / np.sqrt(det)
    try:
        L = logm(mat)
    except Exception:
        return None
    x = float(np.real(L[0,1] + L[1,0])) / 2
    y = float(np.imag(L[1,0] - L[0,1])) / 2
    z = float(np.real(L[0,0] - L[1,1])) / 2
    v = np.array([x, y, z])
    nv = np.linalg.norm(v)
    if nv < 1e-10:
        return None
    return v / nv

# ── Vectorized scan for one manifold ─────────────────────────────
def scan_manifold(rho, words):
    # Build axes
    axes = {}
    for w in words:
        try:
            mat = np.array(rho(w), dtype=complex)
            ax  = matrix_to_axis(mat)
            if ax is not None:
                axes[w] = ax
        except Exception:
            continue
    if len(axes) < 3:
        return None

    valid = list(axes.keys())
    AX    = np.array([axes[w] for w in valid])
    nw    = len(valid)

    # Pairwise angle matrix
    dots  = np.clip(np.abs(AX @ AX.T), 0, 1)
    THETA = np.arccos(dots)

    # All triples
    triple_idx = np.array(list(combinations(range(nw), 3)))
    i, j, k = triple_idx[:,0], triple_idx[:,1], triple_idx[:,2]
    t01 = THETA[i, j]
    t02 = THETA[i, k]
    t12 = THETA[j, k]

    best_fitness = 1e9
    best_rec     = None

    for sigma in SIGMA_VALUES:
        s2  = 2 * sigma**2
        o01 = np.exp(-t01**2 / s2)
        o02 = np.exp(-t02**2 / s2)
        o12 = np.exp(-t12**2 / s2)

        # Proxy filter
        proxy = (np.abs(o01 - PMNS[0,1]) +
                 np.abs(o02 - PMNS[0,2]) +
                 np.abs(o12 - PMNS[1,2]))
        top_idx = np.argpartition(proxy, min(TOP_CANDIDATES, len(proxy)-1))[:TOP_CANDIDATES]

        for ti in top_idx:
            ii, jj, kk = triple_idx[ti]
            O = np.array([[1.0,     o01[ti], o02[ti]],
                          [o01[ti], 1.0,     o12[ti]],
                          [o02[ti], o12[ti], 1.0    ]])
            U, _ = qr(O)
            for a in range(3):
                if U[a,a] < 0: U[:,a] *= -1
            f = float(np.linalg.norm(np.abs(U) - PMNS, 'fro'))
            if f < best_fitness:
                best_fitness = f
                best_rec = {
                    'sigma':   round(float(sigma), 4),
                    'fitness': round(f, 6),
                    'word1':   valid[ii],
                    'word2':   valid[jj],
                    'word3':   valid[kk],
                    'theta12': round(float(np.degrees(THETA[ii,jj])), 2),
                    'theta13': round(float(np.degrees(THETA[ii,kk])), 2),
                    'theta23': round(float(np.degrees(THETA[jj,kk])), 2),
                }

    return best_rec

# ── Main ──────────────────────────────────────────────────────────
def main():
    words = gen_words(MAX_WORDLEN)
    print(f"Word pool: {len(words)}")
    print(f"Scanning manifolds 0-{CENSUS_RANGE-1}")
    print(f"Sigma: [{SIGMA_VALUES[0]:.2f}, {SIGMA_VALUES[-1]:.2f}] in {len(SIGMA_VALUES)} steps")
    print(f"Results -> {RESULTS_FILE}")
    print("=" * 65)

    results   = []
    global_best = 1e9
    start_time  = time.time()

    # Checkpoint resume
    start_idx = 0
    if os.path.exists(RESULTS_FILE):
        try:
            existing = pd.read_csv(RESULTS_FILE)
            if len(existing) > 0:
                start_idx   = int(existing['census_idx'].max()) + 1
                results     = existing.to_dict('records')
                global_best = existing['fitness'].min()
                print(f"Resuming from {start_idx}, best so far: {global_best:.6f}")
        except Exception:
            pass

    for idx in range(start_idx, CENSUS_RANGE):
        t0 = time.time()
        try:
            M   = snappy.OrientableClosedCensus[idx]
            rho = M.polished_holonomy()
            vol = float(M.volume())
            h1  = str(M.homology())
        except Exception as e:
            print(f"[{idx:3d}] ERROR: {e}")
            continue

        rec = scan_manifold(rho, words)
        elapsed = time.time() - t0

        if rec:
            rec.update({
                'census_idx': idx,
                'manifold':   M.name(),
                'volume':     round(vol, 5),
                'homology':   h1,
            })
            results.append(rec)
            if rec['fitness'] < global_best:
                global_best = rec['fitness']

            # Checkpoint
            df = pd.DataFrame(results)
            df.nsmallest(TOP_N, 'fitness').to_csv(RESULTS_FILE, index=False)

        total_elapsed = time.time() - start_time
        done = idx - start_idx + 1
        eta  = (total_elapsed / done) * (CENSUS_RANGE - idx - 1)

        marker = " <<<" if rec and rec['fitness'] < 0.20 else ""
        print(f"[{idx:3d}] {M.name():8s} vol={vol:.4f} h1={h1:6s} "
              f"best={rec['fitness'] if rec else 9.9:.5f} "
              f"global={global_best:.5f} t={elapsed:.1f}s "
              f"ETA={eta/60:.1f}m{marker}", flush=True)

        # Progress file for monitor
        with open(PROGRESS_FILE, 'w') as pf:
            pf.write(f"{idx},{CENSUS_RANGE},{global_best:.6f},{eta:.0f}\n")

    print("=" * 65)
    print("SCAN COMPLETE")
    df = pd.DataFrame(results)
    print(df.nsmallest(10, 'fitness').to_string(index=False))

if __name__ == '__main__':
    main()