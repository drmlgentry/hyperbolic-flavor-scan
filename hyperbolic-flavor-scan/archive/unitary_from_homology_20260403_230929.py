import numpy as np

pi = np.pi

# -----------------------------
# Core math
# -----------------------------
def phases_from_classes(classes, p):
    return [2*pi*k/p for k in classes]

def dft3():
    w = np.exp(2j*np.pi/3)
    F = np.array([
        [1, 1, 1],
        [1, w, w**2],
        [1, w**2, w]
    ], dtype=complex)
    return F / np.sqrt(3)

def build_unitary(classes, p=5):
    phi = phases_from_classes(classes, p)
    D = np.diag(np.exp(1j*np.array(phi)))
    F = dft3()
    return F @ D

def jarlskog(U):
    return np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0]))

def abs_matrix(U):
    return np.abs(U)

# -----------------------------
# Input from your scan
# -----------------------------
m003_classes = [2,4,0]
m006_classes = [2,4,0]

# -----------------------------
# Build matrices
# -----------------------------
U_pmns = build_unitary(m003_classes)
U_ckm  = build_unitary(m006_classes)

# -----------------------------
# Output
# -----------------------------
np.set_printoptions(precision=4, suppress=True)

print("=== PMNS (m003) ===")
print(abs_matrix(U_pmns))
print("J =", jarlskog(U_pmns))

print("\n=== CKM (m006) ===")
print(abs_matrix(U_ckm))
print("J =", jarlskog(U_ckm))

# -----------------------------
# Reference experimental values
# -----------------------------
CKM_exp = np.array([
    [0.974, 0.226, 0.0036],
    [0.226, 0.973, 0.042],
    [0.0087, 0.041, 0.999]
])

PMNS_exp = np.array([
    [0.82, 0.55, 0.15],
    [0.50, 0.58, 0.64],
    [0.28, 0.60, 0.75]
])

print("\n=== Distance to experiment (Frobenius norm) ===")
print("PMNS distance:", np.linalg.norm(abs_matrix(U_pmns) - PMNS_exp))
print("CKM  distance:", np.linalg.norm(abs_matrix(U_ckm)  - CKM_exp))
