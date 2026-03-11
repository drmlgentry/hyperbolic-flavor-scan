import snappy
import numpy as np
import pandas as pd
import itertools
from scipy.linalg import logm, qr

# ----------------------------------------------------------------------
# Copy of matrix_to_axis_vector from refine_flavor_match.py (exact)
# ----------------------------------------------------------------------
def matrix_to_axis_vector(matrix):
    """
    Extract a direction vector from an SL(2,C) matrix.
    Uses the matrix logarithm to get the "axis" direction.
    """
    try:
        # Convert to numpy array with complex entries
        mat = np.array(matrix, dtype=complex)
        
        # Normalize to SL(2,C)
        det = np.linalg.det(mat)
        mat = mat / np.sqrt(det)
        
        # Take matrix logarithm
        log_mat = logm(mat)
        
        # The "direction" is encoded in the off-diagonal structure
        # Extract a characteristic vector
        a = log_mat[0, 0]
        b = log_mat[0, 1]
        c = log_mat[1, 0]
        d = log_mat[1, 1]
        
        # Create a 3-vector from the traceless part
        # Using the Pauli decomposition: M = t*I + x*σ1 + y*σ2 + z*σ3
        x = float(np.real(b + c)) / 2
        y = float(np.imag(c - b)) / 2
        z = float(np.real(a - d)) / 2
        
        vec = np.array([x, y, z])
        norm = np.linalg.norm(vec)
        if norm > 1e-10:
            return vec / norm
        else:
            return np.array([1.0, 0.0, 0.0])
    except:
        return np.array([1.0, 0.0, 0.0])

# ----------------------------------------------------------------------
# Copy of vector_angle from refine_flavor_match.py (uses absolute dot)
# ----------------------------------------------------------------------
def vector_angle(v1, v2):
    """Angle between two unit vectors in radians."""
    dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
    return np.arccos(np.abs(dot))  # Use abs to treat parallel and antiparallel as same

# ----------------------------------------------------------------------
# Build mixing matrix from three words (following evaluate_manifold logic)
# ----------------------------------------------------------------------
def mixing_matrix_from_words(rho, words, sigma, alpha=0.0):
    """
    Returns U (unitary mixing matrix) and its moduli.
    If alpha=0, phases are zero.
    """
    # Compute axes
    axes = [matrix_to_axis_vector(rho(w)) for w in words]
    
    # Compute phases (homology) – simplified: alpha=0 => zero phases
    # (We could compute homology images, but with alpha=0 they don't matter)
    phases = [0.0, 0.0, 0.0]  # all zero
    
    # Build overlap matrix with Gaussian kernel and phases
    mo = np.zeros((3, 3), dtype=complex)
    for r in range(3):
        for c in range(3):
            angle = vector_angle(axes[r], axes[c])
            overlap = np.exp(-(angle**2) / (2 * sigma**2))
            phase_factor = np.exp(1j * (phases[r] - phases[c]))
            mo[r, c] = overlap * phase_factor
    
    # QR orthogonalization
    U, _ = qr(mo)
    return U

# ----------------------------------------------------------------------
# Main scan
# ----------------------------------------------------------------------
# Load m006
M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()

# CKM target (moduli)
CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

SIGMA = 0.49   # best sigma
ALPHA = 0.0    # no phases (focus on moduli)

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

# Collect hyperbolic words
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
        U = mixing_matrix_from_words(rho, list(triple), sigma=SIGMA, alpha=ALPHA)
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
        # Skip problematic triples silently
        pass

    count += 1
    if count % 100 == 0:
        print(f"Processed {count}/{len(triples)} triples...")

# Save to CSV
df = pd.DataFrame(results)
df = df.sort_values('fitness')
df.to_csv('word_triple_results_corrected.csv', index=False)
print("\nWord triple scan complete. Saved to word_triple_results_corrected.csv")
print("Top 10 triples:")
print(df.head(10).to_string())
