import snappy
import numpy as np

# -----------------------------
# Adaptive geodesic extraction
# -----------------------------
def get_lengths(M, n=3, start=0.5, step=0.5, max_L=6.0):
    L = start
    while L <= max_L:
        spec = M.length_spectrum(L)
        lengths = []

        for g in spec:
            try:
                lengths.append(float(g.length.real))
            except:
                continue

        lengths = sorted(set(lengths))

        if len(lengths) >= n:
            return lengths[:n]

        L += step

    return lengths  # fallback (may be < n)

# -----------------------------
# Phase + unitary construction
# -----------------------------
def phases(classes, p):
    return np.array([2*np.pi*k/p for k in classes], dtype=float)

def build_unitary(classes, lengths, p):
    phi = phases(classes, p)
    lengths = np.array(lengths, dtype=float)

    # weights from lengths
    w = np.exp(-lengths/2.0)
    w = w / np.linalg.norm(w)

    U = np.zeros((3,3), dtype=complex)

    for j in range(3):
        for i in range(3):
            U[i,j] = w[j] * np.exp(1j*(phi[j] + 2*np.pi*i*j/3))

    # Orthonormalize
    Q, _ = np.linalg.qr(U)
    return Q

def jarlskog(U):
    return np.imag(
        U[0,0]*np.conj(U[0,1])*
        np.conj(U[1,0])*U[1,1]
    )

# -----------------------------
# Cases
# -----------------------------
cases = [
    ("m003", 1, [2,4,0], 5),
    ("m006", 43, [2,4,0], 5),
]

np.set_printoptions(precision=4, suppress=True)

for name, idx, classes, p in cases:
    print(f"\n=== {name} (idx {idx}) ===")

    M = snappy.OrientableClosedCensus[idx]

    lengths = get_lengths(M)

    print("Lengths:", lengths)

    if len(lengths) < 3:
        print("WARNING: fewer than 3 geodesics found")
        continue

    U = build_unitary(classes, lengths, p)

    print("|U| =")
    print(np.abs(U))

    print("J =", float(jarlskog(U)))
