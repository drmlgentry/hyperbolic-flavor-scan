import snappy
import numpy as np

def get_lengths(M, n=3):
    # Convert SnapPy Number -> float immediately
    spec = M.length_spectrum(0.6)
    lengths = []
    for g in spec:
        try:
            lengths.append(float(g.length.real))
        except:
            continue
        if len(lengths) >= n:
            break
    return lengths

def phases_from_classes(classes, p):
    return [2*np.pi*float(k)/p for k in classes]

def build_unitary(classes, lengths, p):
    phi = phases_from_classes(classes, p)
    lengths = np.array(lengths, dtype=float)

    # Use lengths to define weights (non-uniform)
    w = lengths / np.linalg.norm(lengths)

    # Build columns with distinct structure (avoid rank-1)
    U = np.zeros((3,3), dtype=complex)

    for j in range(3):
        for i in range(3):
            phase = np.exp(1j * (phi[j] + 2*np.pi*i*j/3))
            U[i,j] = w[j] * phase

    # QR orthonormalization
    Q, _ = np.linalg.qr(U)
    return Q

def jarlskog(U):
    return np.imag(
        U[0,0]*np.conj(U[0,1])*
        np.conj(U[1,0])*U[1,1]
    )

cases = [
    ("m003", 1, [2,4,0], 5),
    ("m006", 43, [2,4,0], 5),
]

for name, idx, classes, p in cases:
    print(f"\n=== {name} (idx {idx}) ===")
    M = snappy.OrientableClosedCensus[idx]

    lengths = get_lengths(M)
    print("Lengths:", lengths)

    if len(lengths) < 3:
        print("Not enough geodesics")
        continue

    U = build_unitary(classes, lengths, p)
    print("U =")
    print(np.round(U,4))

    print("J =", float(jarlskog(U)))
