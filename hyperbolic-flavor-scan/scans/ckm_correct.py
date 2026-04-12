import snappy
import numpy as np
from scipy.linalg import logm, qr

def matrix_to_axis_vector(matrix):
    mat = np.array(matrix, dtype=complex)
    det = np.linalg.det(mat)
    mat = mat / np.sqrt(det)
    log_mat = logm(mat)
    a = log_mat[0,0]; b = log_mat[0,1]
    c = log_mat[1,0]; d = log_mat[1,1]
    x = float(np.real(b + c)) / 2
    y = float(np.imag(c - b)) / 2
    z = float(np.real(a - d)) / 2
    vec = np.array([x, y, z])
    norm = np.linalg.norm(vec)
    return vec/norm if norm > 1e-10 else np.array([1.,0.,0.])

def vector_angle(v1, v2):
    dot = np.clip(np.dot(v1,v2), -1., 1.)
    return np.arccos(np.abs(dot))

def mixing_matrix(rho, words, sigma=0.49):
    axes = [matrix_to_axis_vector(rho(w)) for w in words]
    mo = np.zeros((3,3), dtype=complex)
    for r in range(3):
        for c in range(3):
            angle = vector_angle(axes[r], axes[c])
            mo[r,c] = np.exp(-(angle**2)/(2*sigma**2))
    U, _ = qr(mo)
    return U, axes

CKM_PDG = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

PMNS_PDG = np.array([
    [0.8220, 0.5470, 0.1520],
    [0.4300, 0.6110, 0.6640],
    [0.3710, 0.5720, 0.7320]
])

print("=== CKM from m006 ===")
M_ckm = snappy.OrientableClosedCensus[43]
rho_ckm = M_ckm.polished_holonomy()
ckm_words = ['aaB', 'AbA', 'AAb']
U_ckm, axes_ckm = mixing_matrix(rho_ckm, ckm_words)
fit_ckm = np.linalg.norm(np.abs(U_ckm) - CKM_PDG, 'fro')
print(f"Words: {ckm_words}  sigma=0.49")
print(f"Fitness: {fit_ckm:.5f}")
print(f"|U_CKM|:")
for row in np.abs(U_ckm):
    print("  " + "  ".join(f"{x:.5f}" for x in row))
print(f"PDG CKM:")
for row in CKM_PDG:
    print("  " + "  ".join(f"{x:.5f}" for x in row))
print(f"\nAxes (logm/Pauli method):")
for w, n in zip(ckm_words, axes_ckm):
    print(f"  {w}: ({n[0]:.4f}, {n[1]:.4f}, {n[2]:.4f})")

print()
print("=== PMNS from m003 ===")
M_pmns = snappy.OrientableClosedCensus[1]
rho_pmns = M_pmns.polished_holonomy()

# Scan for best PMNS triple
import itertools
letters = ['a','b','A','B']
all_words = [''.join(p) for L in range(1,4)
             for p in itertools.product(letters, repeat=L)]

hyp_words = []
for w in all_words:
    try:
        mat = rho_pmns(w)
        if abs(float(np.abs(np.trace(np.array(mat, dtype=complex))))) > 2.01:
            hyp_words.append(w)
    except:
        pass
print(f"Hyperbolic words: {len(hyp_words)}")

best_fit = 999
best_triple = None
best_U = None
for triple in itertools.combinations(hyp_words, 3):
    try:
        U, _ = mixing_matrix(rho_pmns, list(triple))
        fit = np.linalg.norm(np.abs(U) - PMNS_PDG, 'fro')
        if fit < best_fit:
            best_fit = fit
            best_triple = triple
            best_U = U
    except:
        pass

print(f"Best PMNS triple: {best_triple}")
print(f"Fitness: {best_fit:.5f}")
print(f"|U_PMNS|:")
for row in np.abs(best_U):
    print("  " + "  ".join(f"{x:.5f}" for x in row))
print(f"PDG PMNS:")
for row in PMNS_PDG:
    print("  " + "  ".join(f"{x:.5f}" for x in row))

print()
print("=== Qubit interpretation of logm/Pauli axes ===")
print()
print("The logm/Pauli method extracts from M = exp(log M):")
print("  log M = traceless part => a*sigma_x + b*sigma_y + c*sigma_z")
print("  Axis = (Re(b+c)/2, Im(c-b)/2, Re(a-d)/2) normalized")
print()
print("This is the IMAGINARY part of the matrix logarithm")
print("projected onto su(2) via Pauli basis.")
print()
print("Qubit interpretation:")
print("  This axis IS a Bloch vector -- but from Im(log M), not Re(log M)")
print("  The SU(2) rotation axis of M is exactly this Bloch vector")
print("  The Gaussian kernel sigma=0.49 sets the 'coherence length'")
print("  on the Bloch sphere for the mixing amplitude")
print()
PHI = (1+5**0.5)/2
print(f"sigma = 0.49 rad = {np.degrees(0.49):.2f} deg")
print(f"1/phi^2 = {1/PHI**2:.4f}  vs  sigma = 0.49")
print(f"log(phi) = {np.log(PHI):.4f}  vs  sigma = 0.49")
print(f"sigma * sqrt(2) = {0.49*2**0.5:.4f}")
print(f"sigma / log(phi) = {0.49/np.log(PHI):.4f}")
