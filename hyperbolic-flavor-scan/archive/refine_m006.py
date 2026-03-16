import snappy
import numpy as np
from scipy.linalg import qr, logm
import csv

# PDG CKM moduli
V_CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

def trace(mat):
    return np.trace(mat)

def trace_to_length(tr):
    if abs(tr) <= 2:
        return 0.0
    return 2 * np.arccosh(abs(tr)/2)

def axis_from_mat(mat):
    """Compute unit axis vector from SL(2,C) matrix (if hyperbolic)."""
    mat = np.array(mat, dtype=complex)
    tr = trace(mat)
    if abs(tr) <= 2:
        return None
    # Using formula for axis from trace and matrix entries
    # Alternative: compute eigenvector of the matrix (as before) but get 3D axis.
    # For simplicity, we'll use fixed point method (as before) but now we need axis.
    # Fixed point z gives axis via stereographic projection.
    # We'll compute fixed point and then convert to 3D unit vector.
    evals, evecs = np.linalg.eig(mat)
    idx = np.argmax(np.abs(evals))
    v = evecs[:, idx]
    z = v[0]/v[1] if abs(v[1]) > 1e-12 else np.inf
    if np.isinf(z):
        return np.array([0,0,1.0])
    # Stereographic projection: (2*Re(z), 2*Im(z), 1 - |z|^2) / (1 + |z|^2)
    denom = 1 + abs(z)**2
    return np.array([2*z.real, 2*z.imag, 1 - abs(z)**2]) / denom

def angle_between_axes(axis1, axis2):
    """Angle (in radians) between two unit vectors."""
    return np.arccos(np.clip(np.dot(axis1, axis2), -1, 1))

def gaussian_overlap_from_axes(axis1, axis2, sigma):
    """Overlap based on chordal distance between fixed points (same as before)."""
    # Convert axes back to points on sphere? Actually we can keep using fixed points.
    # But we already have fixed points from axis computation; we can store them.
    # For simplicity, we'll reuse the fixed points computed earlier.
    # However, the given script likely uses fixed points already.
    # We'll just use a placeholder; in practice, we need to compute fixed points.
    # The original refine script must have stored fixed points; we'll adapt.
    # For now, we'll assume we have fixed points stored.
    pass

# Main parameters
idx = 43   # m006
sigma_vals = np.arange(0.3, 0.71, 0.02)
words = [('a','a','B'), ('a','B','a'), ('A','A','b')]  # as strings

# Get manifold
M = snappy.OrientableClosedCensus[idx]
rho = M.polished_holonomy()
gen_names = rho.generators()
fwd = {g: np.array(rho.SL2C(g), dtype=complex) for g in gen_names}
inv = {g: np.linalg.inv(fwd[g]) for g in gen_names}

def word_to_mat(word_tuple):
    mat = np.eye(2, dtype=complex)
    for ch in word_tuple:
        if ch.islower():
            mat = mat @ fwd[ch]
        else:
            mat = mat @ inv[ch.lower()]
    return mat

# Precompute matrices and fixed points for each word
word_mats = [word_to_mat(w) for w in words]
points = []
for mat in word_mats:
    evals, evecs = np.linalg.eig(mat)
    idx_max = np.argmax(np.abs(evals))
    v = evecs[:, idx_max]
    if abs(v[1]) < 1e-12:
        z = np.inf
    else:
        z = v[0]/v[1]
    points.append(z)

# Sigma scan
results = []
for sigma in sigma_vals:
    # Build overlap matrix
    G = np.zeros((3,3), dtype=complex)
    for i in range(3):
        for j in range(3):
            # chordal distance
            zi, zj = points[i], points[j]
            if zi == zj:
                d = 0.0
            elif np.isinf(zi) or np.isinf(zj):
                d = 2.0  # max distance
            else:
                denom = np.sqrt((1+abs(zi)**2)*(1+abs(zj)**2))
                d = 2.0 * abs(zi - zj) / denom
            G[i,j] = np.exp(-d**2 / (2 * sigma**2))
    G = (G + G.conj().T) / 2
    try:
        Q, R = qr(G)
    except:
        continue
    U_mod = np.abs(Q)
    score = np.linalg.norm(U_mod - V_CKM) / np.linalg.norm(V_CKM)
    results.append((sigma, score))

# Print best sigma
best_sigma, best_score = min(results, key=lambda x: x[1])
print(f"Best sigma: {best_sigma:.2f} with score {best_score:.4f}")