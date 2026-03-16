import snappy
import numpy as np
import csv
import sys
from itertools import product

# PDG CKM moduli (central values)
V_CKM_TARGET = np.array([
    [0.97401, 0.22650, 0.00361],
    [0.22650, 0.97359, 0.04153],
    [0.00854, 0.04053, 0.99910]
])

J_TARGET = 3.00e-5

def chordal_dist(z1, z2):
    if z1 == z2:
        return 0.0
    if np.isinf(z1):
        z1 = 1e12
    if np.isinf(z2):
        z2 = 1e12
    denom = np.sqrt((1 + abs(z1)**2) * (1 + abs(z2)**2))
    return 2.0 * abs(z1 - z2) / denom

def gaussian_overlap(z1, z2, sigma=0.5):
    d = chordal_dist(z1, z2)
    return np.exp(-d**2 / (2 * sigma**2))

def fixed_point(mat):
    mat = np.array(mat, dtype=complex)
    evals, evecs = np.linalg.eig(mat)
    idx = np.argmax(np.abs(evals))
    v = evecs[:, idx]
    if abs(v[1]) < 1e-12:
        return np.inf
    return v[0] / v[1]

def jarlskog(U):
    return np.imag(U[0,0] * U[1,1] * np.conj(U[0,1]) * np.conj(U[1,0]))

def build_mixing_matrix(points, sigma):
    n = len(points)
    G = np.array([[gaussian_overlap(points[i], points[j], sigma) for j in range(n)]
                  for i in range(n)], dtype=complex)
    G = (G + G.conj().T) / 2
    try:
        Q, R = np.linalg.qr(G)
        return Q
    except:
        return None

def fitness(U_mod, J):
    mod_error = np.linalg.norm(U_mod - V_CKM_TARGET) / np.linalg.norm(V_CKM_TARGET)
    j_error = abs(J - J_TARGET) / J_TARGET if J_TARGET != 0 else abs(J)
    return mod_error + 0.1 * j_error

def trace(mat):
    return np.trace(mat)

def trace_to_length(tr):
    """Return translation length if hyperbolic, else 0."""
    tr = complex(tr)
    if abs(tr) <= 2:
        return 0.0
    return 2 * np.arccosh(abs(tr) / 2)

def generate_words(num_gens, max_len):
    """Generate all non‑empty words up to max_len using forward and inverse letters.
       Returns list of tuples of integers: positive for forward, negative for inverse.
    """
    letters = list(range(1, num_gens+1)) + list(range(-num_gens, 0))
    words = []
    for length in range(1, max_len+1):
        for prod in product(letters, repeat=length):
            words.append(prod)
    return words

def word_to_matrix(word, fwd_mats, inv_mats):
    """Return matrix product (right‑multiplying) for a word.
       fwd_mats[i] corresponds to generator i+1.
       inv_mats[i] is the inverse of fwd_mats[i].
    """
    mat = np.eye(2, dtype=complex)
    for w in word:
        if w > 0:
            mat = mat @ fwd_mats[w-1]
        else:
            mat = mat @ inv_mats[-w-1]
    return mat

def get_three_shortest_matrices(M, max_len=4):
    """Return three numpy arrays (SL(2,C) matrices) corresponding to the
       three shortest hyperbolic elements in the fundamental group."""
    try:
        rho = M.polished_holonomy()
        gen_names = rho.generators()
        if len(gen_names) < 2:
            print(f"  Too few generators: {len(gen_names)}")
            return None
        # Use the first two generators (most manifolds have exactly two)
        fwd_mats = [np.array(rho.SL2C(gen_names[0]), dtype=complex),
                    np.array(rho.SL2C(gen_names[1]), dtype=complex)]
        inv_mats = [np.linalg.inv(m) for m in fwd_mats]

        words = generate_words(2, max_len)
        traces = []
        for w in words:
            try:
                mat = word_to_matrix(w, fwd_mats, inv_mats)
                tr = trace(mat)
                abs_tr = abs(tr)
                if abs_tr > 2:   # hyperbolic
                    traces.append((abs_tr, w, mat))
            except:
                continue

        if len(traces) < 3:
            print(f"  Only {len(traces)} hyperbolic words found (need 3).")
            return None

        # Sort by absolute trace
        traces.sort(key=lambda x: x[0])

        # Deduplicate by absolute trace (rounded)
        unique = []
        seen = set()
        for abs_tr, w, mat in traces:
            rounded = round(abs_tr, 8)
            if rounded not in seen:
                seen.add(rounded)
                unique.append(mat)
            if len(unique) >= 3:
                break

        if len(unique) < 3:
            print(f"  Insufficient distinct hyperbolic traces.")
            return None

        return unique[:3]

    except Exception as e:
        print(f"  Error in get_three_shortest_matrices: {e}")
        return None

def scan_manifolds(limit=50, sigma=0.5, max_word_len=4, outfile='scan_results.csv'):
    results = []
    with open(outfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'name', 'volume', 'fitness', 'jarlskog',
                         'u11', 'u12', 'u13', 'u21', 'u22', 'u23', 'u31', 'u32', 'u33'])

    for i in range(limit):
        M = snappy.OrientableClosedCensus[i]
        name = M.name()
        try:
            vol = float(M.volume())
        except:
            vol = -1.0
        print(f"Processing {name} (vol={vol:.4f})...")

        try:
            matrices = get_three_shortest_matrices(M, max_len=max_word_len)
            if matrices is None or len(matrices) < 3:
                print("  Could not get three shortest matrices, skipping.")
                continue

            points = []
            for g in matrices:
                try:
                    pt = fixed_point(g)
                    points.append(pt)
                except:
                    points.append(np.inf)

            U = build_mixing_matrix(points, sigma)
            if U is None:
                print("  Failed to build mixing matrix.")
                continue

            U_mod = np.abs(U)
            J = jarlskog(U)
            score = fitness(U_mod, J)

            results.append({
                'name': name,
                'volume': vol,
                'score': score,
                'J': J,
                'U_mod': U_mod
            })

            with open(outfile, 'a', newline='') as f:
                writer = csv.writer(f)
                flat = [i, name, f"{vol:.6f}", f"{score:.6f}", f"{J:.2e}"]
                flat.extend([f"{x:.4f}" for row in U_mod for x in row])
                writer.writerow(flat)

            print(f"  score={score:.4f}, J={J:.2e}")
        except Exception as e:
            print(f"  Error: {e}")
            continue

    results.sort(key=lambda x: x['score'])
    return results

if __name__ == "__main__":
    sigma = 0.5
    max_word_len = 4
    outfile = 'scan_results_shortest.csv'
    if len(sys.argv) > 1:
        try:
            sigma = float(sys.argv[1])
        except:
            print(f"Invalid sigma value '{sys.argv[1]}', using default {sigma}")
    if len(sys.argv) > 2:
        try:
            max_word_len = int(sys.argv[2])
        except:
            print(f"Invalid max_word_len '{sys.argv[2]}', using default {max_word_len}")
    if len(sys.argv) > 3:
        outfile = sys.argv[3]

    print(f"Scanning first 50 closed manifolds with sigma={sigma}, max_word_len={max_word_len}...\n")
    results = scan_manifolds(limit=50, sigma=sigma, max_word_len=max_word_len, outfile=outfile)
    print("\nTop 5 candidates:")
    for r in results[:5]:
        print(f"{r['name']:20s} vol {r['volume']:.4f}  score {r['score']:.4f}  J {r['J']:.2e}")