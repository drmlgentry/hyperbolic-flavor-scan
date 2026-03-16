import snappy
import numpy as np
import pandas as pd
import itertools
from scipy.linalg import logm, qr

# ----------------------------------------------------------------------
# Helper: axis from SL(2,C) matrix (converts PARI -> NumPy complex)
# ----------------------------------------------------------------------
def matrix_to_axis_vector(mat):
    try:
        mat_np = np.array(mat, dtype=complex)
        det = np.linalg.det(mat_np)
        mat_np = mat_np / np.sqrt(det)
        L = logm(mat_np)
        n = np.array([L[1,0].imag, L[0,1].imag, L[0,0].imag])
        norm = np.linalg.norm(n)
        if norm < 1e-12:
            n = np.array([L[1,0].real, L[0,1].real, L[0,0].real])
            norm = np.linalg.norm(n)
        return n / norm
    except:
        return np.array([1.0, 0.0, 0.0])

# ----------------------------------------------------------------------
# Helper: compute mixing matrix from three words
# ----------------------------------------------------------------------
def mixing_matrix_from_words(rho, words, sigma):
    axes = []
    for w in words:
        mat = rho(w)
        axis = matrix_to_axis_vector(mat)
        axes.append(axis)

    angles = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            cos = np.dot(axes[i], axes[j])
            cos = np.clip(cos, -1.0, 1.0)
            angles[i,j] = np.arccos(cos)

    K = np.exp(- (angles**2) / (2 * sigma**2))
    Q, _ = qr(K)
    return Q, axes, angles

# ----------------------------------------------------------------------
# Main: word triple scan for m006 – all hyperbolic words length ≤ 3
# ----------------------------------------------------------------------
M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()

CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

SIGMA = 0.49   # fixed from best result

# Generate all words up to length 3
letters = ['a', 'b', 'A', 'B']
all_words = []
for L in range(1, 4):
    for prod in itertools.product(letters, repeat=L):
        all_words.append(''.join(prod))
print(f"Generated {len(all_words)} words up to length 3")

def is_hyperbolic(mat):
    try:
        tr = float(np.abs(np.trace(mat)))
        return tr > 2.01
    except:
        return False

# Collect hyperbolic words (length ≤ 3)
hyperbolic_words = []
for w in all_words:
    try:
        mat = rho(w)
        if is_hyperbolic(mat):
            hyperbolic_words.append(w)
    except Exception:
        continue

print(f"Found {len(hyperbolic_words)} hyperbolic words of length ≤ 3")

# Generate all triples
triples = list(itertools.combinations(hyperbolic_words, 3))
print(f"Testing {len(triples)} triples")

results = []
count = 0
for triple in triples:
    try:
        U, axes, angles = mixing_matrix_from_words(rho, list(triple), sigma=SIGMA)
        diff = np.abs(U) - CKM
        fitness = np.linalg.norm(diff, 'fro')
        results.append({
            'words': ' '.join(triple),
            'fitness': fitness,
            'U00': np.abs(U[0,0]), 'U01': np.abs(U[0,1]), 'U02': np.abs(U[0,2]),
            'U10': np.abs(U[1,0]), 'U11': np.abs(U[1,1]), 'U12': np.abs(U[1,2]),
            'U20': np.abs(U[2,0]), 'U21': np.abs(U[2,1]), 'U22': np.abs(U[2,2])
        })
    except Exception as e:
        # Silently skip problematic triples
        pass

    count += 1
    if count % 1000 == 0:
        print(f"Processed {count}/{len(triples)} triples...")

# Save to CSV
df = pd.DataFrame(results)
df = df.sort_values('fitness')
df.to_csv('word_triple_results_full.csv', index=False)
print("\nWord triple scan complete. Saved to word_triple_results_full.csv")
print("Top 10 triples:")
print(df.head(10).to_string())
