import snappy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
# Main sigma sensitivity scan
# ----------------------------------------------------------------------
# Load manifold m006 (index 43)
M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()

# Fixed words (best triple)
words = ['aaB', 'aBa', 'AAb']

# CKM target (moduli)
CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

# Sigma range
sigmas = np.arange(0.48, 0.501, 0.001)

results = []
for sigma in sigmas:
    U, axes, angles = mixing_matrix_from_words(rho, words, sigma)
    diff = np.abs(U) - CKM
    fitness = np.linalg.norm(diff, 'fro')
    row = [sigma] + np.abs(U).flatten().tolist() + [fitness]
    results.append(row)

# Save to CSV
columns = ['sigma'] + [f'U{i}{j}' for i in range(3) for j in range(3)] + ['fitness']
df = pd.DataFrame(results, columns=columns)
df.to_csv('sigma_sensitivity.csv', index=False)
print("Sigma sensitivity scan complete. Saved to sigma_sensitivity.csv")

# Plot each entry vs sigma
plt.figure(figsize=(10,6))
for i in range(3):
    for j in range(3):
        plt.plot(df['sigma'], df[f'U{i}{j}'], label=f'|U_{i}{j}|')
plt.xlabel('σ')
plt.ylabel('Matrix entry')
plt.legend()
plt.title('Sigma Sensitivity of Mixing Matrix Entries (m006, words = aaB,aBa,AAb)')
plt.grid(True)
plt.savefig('sigma_sensitivity.png')
plt.show()

# Print max variation
print("\nMaximum variation over sigma range:")
for col in [f'U{i}{j}' for i in range(3) for j in range(3)]:
    variation = df[col].max() - df[col].min()
    print(f"{col}: {variation:.6f}")
