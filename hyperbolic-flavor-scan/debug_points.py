import snappy
import numpy as np

def chordal_dist(z1, z2):
    if z1 == z2:
        return 0.0
    if np.isinf(z1):
        z1 = 1e12
    if np.isinf(z2):
        z2 = 1e12
    denom = np.sqrt((1 + abs(z1)**2) * (1 + abs(z2)**2))
    return 2.0 * abs(z1 - z2) / denom

def gaussian_overlap(z1, z2, sigma=0.5):
    d = chordal_dist(z1, z2)
    return np.exp(-d**2 / (2 * sigma**2))

def fixed_point(mat):
    # Ensure matrix is a proper complex numpy array
    mat = np.array(mat, dtype=complex)
    evals, evecs = np.linalg.eig(mat)
    idx = np.argmax(np.abs(evals))
    v = evecs[:, idx]
    if abs(v[1]) < 1e-12:
        return np.inf
    return v[0] / v[1]

def build_mixing_matrix(points, sigma):
    n = len(points)
    G = np.array([[gaussian_overlap(points[i], points[j], sigma) for j in range(n)]
                  for i in range(n)], dtype=complex)
    G = (G + G.conj().T) / 2
    try:
        Q, R = np.linalg.qr(G)
        return Q
    except:
        return None

def get_matrices(M):
    try:
        rho = M.polished_holonomy()
        gen_names = rho.generators()
        if len(gen_names) < 2:
            return None
        # Convert SnapPy matrices to numpy arrays with complex dtype
        m1 = np.array(rho.SL2C(gen_names[0]), dtype=complex)
        m2 = np.array(rho.SL2C(gen_names[1]), dtype=complex)
        matrices = [m1, m2]
        if len(gen_names) < 3:
            m3 = m1 @ m2
            matrices.append(m3)
        else:
            m3 = np.array(rho.SL2C(gen_names[2]), dtype=complex)
            matrices.append(m3)
        return matrices
    except Exception as e:
        print(f"  Error in get_matrices: {e}")
        return None

print("Index  Manifold  Fixed Points (z1, z2, z3)  |U_mod[0,0]  |U_mod[0,1]")
print("-" * 80)
for i in range(5):
    M = snappy.OrientableClosedCensus[i]
    name = M.name()
    mats = get_matrices(M)
    if mats is None:
        print(f"{i:3d}  {name:8s}  - no matrices -")
        continue
    points = []
    for m in mats:
        try:
            pt = fixed_point(m)
            points.append(pt)
        except Exception as e:
            print(f"    fixed_point error: {e}")
            points.append(np.inf)
    U = build_mixing_matrix(points, sigma=0.5)
    if U is None:
        print(f"{i:3d}  {name:8s}  - mixing matrix failed -")
        continue
    U_mod = np.abs(U)
    # Print points with limited precision
    p_str = " ".join([f"{p:8.3f}" if not np.isinf(p) else "inf" for p in points])
    print(f"{i:3d}  {name:8s}  {p_str}  {U_mod[0,0]:.4f}       {U_mod[0,1]:.4f}")