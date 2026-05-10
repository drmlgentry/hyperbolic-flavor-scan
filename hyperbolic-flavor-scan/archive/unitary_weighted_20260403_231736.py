import numpy as np

pi = np.pi

# -----------------------------
# Inputs (from your scan)
# -----------------------------
m003_classes = [2,4,0]
m006_classes = [2,4,0]

# You MUST replace these with real geodesic lengths later
# For now: use distinct vs degenerate patterns
m003_lengths = [1.0, 1.2, 1.4]   # non-degenerate
m006_lengths = [1.0, 1.0, 1.0]   # degenerate

# -----------------------------
# Core functions
# -----------------------------
def phases(classes, p):
    return np.array([2*pi*k/p for k in classes])

def build_unitary(classes, lengths, p=5):
    phi = phases(classes, p)
    w = np.exp(-np.array(lengths)/2.0)

    # build columns explicitly
    U = np.zeros((3,3), dtype=complex)
    for j in range(3):
        col = w[j] * np.exp(1j*phi[j]) * np.ones(3)
        col = col / np.linalg.norm(col)
        U[:,j] = col

    return U

def jarlskog(U):
    return np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0]))

def absU(U):
    return np.abs(U)

# -----------------------------
# Build matrices
# -----------------------------
U_pmns = build_unitary(m003_classes, m003_lengths)
U_ckm  = build_unitary(m006_classes, m006_lengths)

np.set_printoptions(precision=4, suppress=True)

print("=== PMNS (weighted) ===")
print(absU(U_pmns))
print("J =", jarlskog(U_pmns))

print("\n=== CKM (weighted) ===")
print(absU(U_ckm))
print("J =", jarlskog(U_ckm))
