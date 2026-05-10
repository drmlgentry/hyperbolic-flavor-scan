import snappy
import numpy as np

pi = np.pi

# -----------------------------
# INPUT: your word triples
# -----------------------------
cases = [
    ("m003", 1, ["A","B","aaB"], [2,4,0]),
    ("m006", 43, ["a","B","aab"], [2,4,0]),
]

# -----------------------------
# helpers
# -----------------------------
def phases(classes, p):
    return np.array([2*pi*k/p for k in classes])

def get_geodesics(M, cutoff=3.5):
    geos = M.length_spectrum(cutoff=cutoff)
    # return real lengths only
    return [g.length.real for g in geos]

def assign_lengths(word_list, lengths):
    """
    crude matching:
    assign shortest distinct geodesics to words
    """
    lengths = sorted(lengths)
    if len(lengths) < len(word_list):
        return None
    return lengths[:len(word_list)]

def build_unitary(classes, lengths, p=5):
    phi = phases(classes, p)
    w = np.exp(-np.array(lengths)/2.0)

    U = np.zeros((3,3), dtype=complex)

    for j in range(3):
        col = w[j] * np.exp(1j*phi[j]) * np.ones(3)
        col /= np.linalg.norm(col)
        U[:,j] = col

    return U

def jarlskog(U):
    return np.imag(U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0]))

def absU(U):
    return np.abs(U)

# -----------------------------
# main loop
# -----------------------------
np.set_printoptions(precision=4, suppress=True)

for name, idx, words, classes in cases:
    print(f"\n=== {name} (idx {idx}) ===")

    M = snappy.OrientableClosedCensus[idx]

    # get geodesics
    lengths = get_geodesics(M)

    if not lengths:
        print("No geodesics found")
        continue

    # assign lengths to words
    assigned = assign_lengths(words, lengths)

    if assigned is None:
        print("Not enough geodesics")
        continue

    print("Words:", words)
    print("Assigned lengths:", assigned)

    # build matrix
    U = build_unitary(classes, assigned)

    print("|U| =")
    print(absU(U))
    print("J =", jarlskog(U))
