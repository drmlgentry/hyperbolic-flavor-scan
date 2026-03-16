import snappy
import numpy as np
import pandas as pd
import itertools
from scipy.linalg import logm, qr

# ----------------------------------------------------------------------
# Helper: axis from SL(2,C) matrix (converts PARI -> NumPy complex)
# ----------------------------------------------------------------------
def matrix_to_axis_vector(mat):
    """
    Extract a unit direction vector from an SL(2,C) matrix.
    Uses the matrix logarithm and Pauli decomposition.
    """
    try:
        # Convert to NumPy complex array (essential!)
        mat_np = np.array(mat, dtype=complex)

        # Ensure determinant is 1 (normalize if needed)
        det = np.linalg.det(mat_np)
        mat_np = mat_np / np.sqrt(det)

        # Matrix logarithm
        L = logm(mat_np)

        # Pauli decomposition: L = (i/2)(x σx + y σy + z σz) + (η/2) I
        # The axis direction is proportional to ( Im(L[1,0]), Im(L[0,1]), Im(L[0,0]) )
        n = np.array([L[1,0].imag, L[0,1].imag, L[0,0].imag])
        norm = np.linalg.norm(n)
        if norm < 1e-12:
            # Fallback (should not happen for hyperbolic elements)
            n = np.array([L[1,0].real, L[0,1].real, L[0,0].real])
            norm = np.linalg.norm(n)
        return n / norm
    except Exception as e:
        print(f"Warning: axis extraction failed: {e}")
        return np.array([1.0, 0.0, 0.0])

# ----------------------------------------------------------------------
# Helper: compute mixing matrix from three words
# ----------------------------------------------------------------------
def mixing_matrix_from_words(rho, words, sigma):
    """
    rho : holonomy function (word -> SL(2,C) matrix)
    words : list of three strings
    sigma : Gaussian width
    Returns: U (3x3 unitary), axes (list of vectors), angles (3x3 angle matrix)
    """
    axes = []
    for w in words:
        mat = rho(w)
        axis = matrix_to_axis_vector(mat)
        axes.append(axis)

    # Angle matrix
    angles = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            cos = np.dot(axes[i], axes[j])
            cos = np.clip(cos, -1.0, 1.0)
            angles[i,j] = np.arccos(cos)

    # Gaussian overlap
    K = np.exp(- (angles**2) / (2 * sigma**2))

    # QR orthogonalization -> unitary mixing matrix
    Q, R = qr(K)
    return Q, axes, angles

# ----------------------------------------------------------------------
# Main: word triple scan for m006
# ----------------------------------------------------------------------
# Load manifold m006 (index 43)
M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()

# CKM target (moduli)
CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

# Fixed sigma from best result
SIGMA = 0.49

# Generate all words up to length 4 over a,b,A,B
letters = ['a', 'b', 'A', 'B']
all_words = []
for L in range(1, 5):
    for prod in itertools.product(letters, repeat=L):
        all_words.append(''.join(prod))
print(f"Generated {len(all_words)} words up to length 4")

def is_hyperbolic(mat):
    """Return True if matrix corresponds to a hyperbolic (loxodromic) element."""
    try:
        tr = float(np.abs(np.trace(mat)))
        # Rough criterion: |trace| > 2 indicates hyperbolic
        return tr > 2.01
    except:
        return False

# Find all hyperbolic words
hyperbolic_words = []
for w in all_words:
    try:
        mat = rho(w)
        if is_hyperbolic(mat):
            hyperbolic_words.append(w)
    except Exception:
        continue

print(f"Found {len(hyperbolic_words)} hyperbolic words up to length 4")

# Sort by length (shortest first) and take top 10
hyperbolic_words.sort(key=len)
candidate_words = hyperbolic_words[:10]
print(f"Top 10 shortest hyperbolic words: {candidate_words}")

# Generate all triples
triples = list(itertools.combinations(candidate_words, 3))
print(f"Testing {len(triples)} triples")

results = []
for triple in triples:
    try:
        U, axes, angles = mixing_matrix_from_words(rho, list(triple), sigma=SIGMA)
        diff = np.abs(U) - CKM
        fitness = np.linalg.norm(diff, 'fro')
        row = {
            'words': ' '.join(triple),
            'fitness': fitness,
            'U00': np.abs(U[0,0]), 'U01': np.abs(U[0,1]), 'U02': np.abs(U[0,2]),
            'U10': np.abs(U[1,0]), 'U11': np.abs(U[1,1]), 'U12': np.abs(U[1,2]),
            'U20': np.abs(U[2,0]), 'U21': np.abs(U[2,1]), 'U22': np.abs(U[2,2])
        }
        results.append(row)
    except Exception as e:
        print(f"Error with triple {triple}: {e}")
        continue

# Save to CSV
df = pd.DataFrame(results)
df = df.sort_values('fitness')
df.to_csv('word_triple_results.csv', index=False)
print("\nWord triple scan complete. Saved to word_triple_results.csv")
print("Top 5 triples:")
print(df.head(5).to_string())
