import snappy
import numpy as np
from scipy.linalg import logm, qr
from itertools import combinations

def matrix_to_axis_vector(matrix):
    mat = np.array(matrix, dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    x = float(np.real(L[0,1]+L[1,0]))/2
    y = float(np.imag(L[1,0]-L[0,1]))/2
    z = float(np.real(L[0,0]-L[1,1]))/2
    v = np.array([x,y,z])
    n = np.linalg.norm(v)
    return v/n if n>1e-10 else np.array([1.,0.,0.])

def vector_angle(v1,v2):
    return np.arccos(np.clip(np.abs(np.dot(v1,v2)),-1.,1.))

def scan_best_triple(rho, sigma, target, max_len=3):
    import itertools
    letters = ['a','b','A','B']
    all_words = [''.join(p) for L in range(1,max_len+1)
                 for p in itertools.product(letters,repeat=L)]
    hyp = []
    for w in all_words:
        try:
            mat = np.array(rho(w), dtype=complex)
            if abs(float(np.abs(np.trace(mat)))) > 2.01:
                hyp.append(w)
        except: pass

    best_fit = 999
    best_triple = None
    best_U = None
    for triple in itertools.combinations(hyp, 3):
        try:
            axes = [matrix_to_axis_vector(rho(w)) for w in triple]
            mo = np.zeros((3,3),dtype=complex)
            for i in range(3):
                for j in range(3):
                    mo[i,j] = np.exp(-vector_angle(axes[i],axes[j])**2/(2*sigma**2))
            U,_ = qr(mo)
            fit = np.linalg.norm(np.abs(U)-target,'fro')
            if fit < best_fit:
                best_fit = fit
                best_triple = triple
                best_U = U
        except: pass
    return best_triple, best_fit, best_U

CKM = np.array([[0.97427,0.22536,0.00355],
                [0.22522,0.97339,0.04108],
                [0.00886,0.04050,0.99914]])

PMNS = np.array([[0.822,0.547,0.152],
                 [0.430,0.611,0.664],
                 [0.371,0.572,0.732]])

print("=== Covering tower: best CKM/PMNS fits ===")
print()

manifolds = [
    (43,  "m006 base (CKM)",  CKM,  0.488),
    (1,   "m003 base (PMNS)", PMNS, 0.500),
    (39,  "m003 deg-2 cover", PMNS, 0.500),
    (238, "m003 deg-3 cover", PMNS, 0.500),
]

for idx, label, target, sigma in manifolds:
    print(f"--- {label} (index {idx}) ---")
    M = snappy.OrientableClosedCensus[idx]
    print(f"vol={float(M.volume()):.4f} H1={M.homology()}")
    rho = M.polished_holonomy()
    triple, fit, U = scan_best_triple(rho, sigma, target)
    print(f"Best triple: {triple}")
    print(f"Fitness: {fit:.5f}")
    if U is not None:
        print(f"|U|:")
        for row in np.abs(U):
            print("  " + "  ".join(f"{x:.4f}" for x in row))
    print()
