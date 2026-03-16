import snappy
import numpy as np
import csv
import sys
import re
from itertools import product
from scipy.linalg import qr

# PDG CKM moduli (central values)
V_CKM_TARGET = np.array([
    [0.97401, 0.22650, 0.00361],
    [0.22650, 0.97359, 0.04153],
    [0.00854, 0.04053, 0.99910]
])

J_TARGET = 3.00e-5

# ---------- Geometry helpers ----------
def chordal_dist(z1, z2):
    if z1 == z2:
        return 0.0
    if np.isinf(z1):
        z1 = 1e12
    if np.isinf(z2):
        z2 = 1e12
    denom = np.sqrt((1 + abs(z1)**2) * (1 + abs(z2)**2))
    return 2.0 * abs(z1 - z2) / denom

def gaussian_amplitude(z1, z2, sigma):
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

def trace(mat):
    return np.trace(mat)

# ---------- Word enumeration ----------
def generate_words(num_gens, max_len):
    letters = list(range(1, num_gens+1)) + list(range(-num_gens, 0))
    words = []
    for L in range(1, max_len+1):
        for prod in product(letters, repeat=L):
            words.append(prod)
    return words

def word_to_matrix(word, fwd_mats, inv_mats):
    mat = np.eye(2, dtype=complex)
    for w in word:
        if w > 0:
            mat = mat @ fwd_mats[w-1]
        else:
            mat = mat @ inv_mats[-w-1]
    return mat

# ---------- Homology and U(1) phases ----------
def word_to_homology_vector(word, gen_maps):
    """
    Convert a word (tuple of integers) to a vector in H1(M, Z) (free part).
    gen_maps: list of integer vectors for each generator (from abelianization).
    If free rank = 0, return empty array.
    """
    if len(gen_maps[0]) == 0:
        return np.array([])
    vec = np.zeros(len(gen_maps[0]), dtype=int)
    for w in word:
        if w > 0:
            vec += gen_maps[w-1]
        else:
            vec -= gen_maps[-w-1]
    return vec

def assign_u1_phases(words, gen_maps, free_angles):
    """Return list of phases (complex numbers) for each word."""
    phases = []
    for w in words:
        vec = word_to_homology_vector(w, gen_maps)
        if vec.size == 0:
            phases.append(1.0)
        else:
            angle = np.dot(vec, free_angles)
            phases.append(np.exp(1j * angle))
    return phases

# ---------- Complex overlap and mixing ----------
def complex_overlap(z1, z2, phase1, phase2, sigma):
    amp = gaussian_amplitude(z1, z2, sigma)
    return amp * phase1 * np.conj(phase2)

def build_complex_mixing_matrix(points, phases, sigma):
    n = len(points)
    G = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            G[i,j] = complex_overlap(points[i], points[j], phases[i], phases[j], sigma)
    G = (G + G.conj().T) / 2
    try:
        Q, R = qr(G)
        return Q
    except:
        return None

def jarlskog(U):
    return np.imag(U[0,0] * U[1,1] * np.conj(U[0,1]) * np.conj(U[1,0]))

def fitness(U_mod, J):
    mod_error = np.linalg.norm(U_mod - V_CKM_TARGET) / np.linalg.norm(V_CKM_TARGET)
    j_error = abs(J - J_TARGET) / J_TARGET if J_TARGET != 0 else abs(J)
    return mod_error + 0.1 * j_error

# ---------- Main scanning function ----------
def scan_manifolds_cp(limit=50, sigma=0.5, max_word_len=4, u1_scale=0.01, outfile='scan_cp.csv'):
    results = []
    with open(outfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['index','name','volume','fitness','jarlskog',
                         'u11','u12','u13','u21','u22','u23','u31','u32','u33'])

    for i in range(limit):
        M = snappy.OrientableClosedCensus[i]
        name = M.name()
        try:
            vol = float(M.volume())
        except:
            vol = -1.0
        print(f"Processing {name} (vol={vol:.4f})...")

        try:
            # 1. Get SL(2,C) matrices
            rho = M.polished_holonomy()
            gen_names = rho.generators()
            if len(gen_names) < 2:
                print("  Too few generators.")
                continue
            fwd_mats = [np.array(rho.SL2C(gen_names[0]), dtype=complex),
                        np.array(rho.SL2C(gen_names[1]), dtype=complex)]
            inv_mats = [np.linalg.inv(m) for m in fwd_mats]

            # 2. Enumerate words and find three shortest hyperbolic
            words_list = generate_words(2, max_word_len)
            candidates = []
            for w in words_list:
                try:
                    mat = word_to_matrix(w, fwd_mats, inv_mats)
                    tr = trace(mat)
                    if abs(tr) > 2:
                        candidates.append((abs(tr), w, mat))
                except:
                    continue
            candidates.sort(key=lambda x: x[0])
            # Deduplicate by rounded absolute trace
            unique_words = []
            seen_tr = set()
            for abs_tr, w, mat in candidates:
                rounded = round(abs_tr, 8)
                if rounded not in seen_tr:
                    seen_tr.add(rounded)
                    unique_words.append(w)
                if len(unique_words) >= 3:
                    break
            if len(unique_words) < 3:
                print("  Insufficient distinct hyperbolic words.")
                continue

            # 3. Compute fixed points for the three words
            points = []
            for w in unique_words:
                mat = word_to_matrix(w, fwd_mats, inv_mats)
                try:
                    pt = fixed_point(mat)
                except:
                    pt = np.inf
                points.append(pt)

            # 4. Homology and U(1) phases (with fallback)
            try:
                G = M.fundamental_group()
                abel = G.abelianization()   # list of integer vectors for each generator
                if abel and len(abel[0]) > 0:
                    free_rank = len(abel[0])
                    np.random.seed(int(vol*1000) if vol>0 else 0)
                    free_angles = np.random.uniform(-u1_scale, u1_scale, free_rank)
                    phases = assign_u1_phases(unique_words, abel, free_angles)
                else:
                    phases = [1.0, 1.0, 1.0]
            except Exception as e:
                print(f"  Homology/abelianization failed: {e}, using real kernel")
                phases = [1.0, 1.0, 1.0]

            # 5. Build complex mixing matrix
            U = build_complex_mixing_matrix(points, phases, sigma)
            if U is None:
                print("  Failed to build mixing matrix.")
                continue

            U_mod = np.abs(U)
            J = jarlskog(U)
            score = fitness(U_mod, J)

            results.append({'name':name, 'volume':vol, 'score':score, 'J':J, 'U_mod':U_mod})

            # Save to CSV
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
    max_len = 4
    outfile = 'scan_cp.csv'
    u1_scale = 0.01
    if len(sys.argv) > 1:
        sigma = float(sys.argv[1])
    if len(sys.argv) > 2:
        max_len = int(sys.argv[2])
    if len(sys.argv) > 3:
        outfile = sys.argv[3]
    if len(sys.argv) > 4:
        u1_scale = float(sys.argv[4])
    print(f"Scanning with σ={sigma}, max_word_len={max_len}, U(1) scale={u1_scale}...\n")
    results = scan_manifolds_cp(limit=50, sigma=sigma, max_word_len=max_len,
                                u1_scale=u1_scale, outfile=outfile)
    print("\nTop 5 candidates:")
    for r in results[:5]:
        print(f"{r['name']:20s} vol {r['volume']:.4f}  score {r['score']:.4f}  J {r['J']:.2e}")