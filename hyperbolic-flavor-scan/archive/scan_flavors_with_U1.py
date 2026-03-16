import snappy
import numpy as np
import csv
import sys
import random

# PDG CKM moduli (central values)
V_CKM_TARGET = np.array([
    [0.97401, 0.22650, 0.00361],
    [0.22650, 0.97359, 0.04153],
    [0.00854, 0.04053, 0.99910]
])

J_TARGET = 3.00e-5

# Scale for U(1) phases (radians) – adjust to get realistic Jarlskog
U1_ANGLE_SCALE = 0.01   # small, to mimic observed CP violation

def chordal_dist(z1, z2):
    if z1 == z2:
        return 0.0
    if np.isinf(z1):
        z1 = 1e12
    if np.isinf(z2):
        z2 = 1e12
    denom = np.sqrt((1 + abs(z1)**2) * (1 + abs(z2)**2))
    return 2.0 * abs(z1 - z2) / denom

def gaussian_amplitude(z1, z2, sigma=0.5):
    d = chordal_dist(z1, z2)
    return np.exp(-d**2 / (2 * sigma**2))

def fixed_point(mat):
    mat = np.array(mat, dtype=complex)
    evals, evecs = np.linalg.eig(mat)
    idx = np.argmax(np.abs(evals))
    v = evecs[:, idx]
    if abs(v[1]) < 1e-12:
        return np.inf
    return v[0] / v[1]

def jarlskog(U):
    return np.imag(U[0,0] * U[1,1] * np.conj(U[0,1]) * np.conj(U[1,0]))

def build_mixing_matrix(points, phases, sigma):
    n = len(points)
    G = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            amp = gaussian_amplitude(points[i], points[j], sigma)
            phase = phases[i] * np.conj(phases[j])
            G[i, j] = amp * phase
    G = (G + G.conj().T) / 2
    try:
        Q, R = np.linalg.qr(G)
        return Q
    except:
        return None

def fitness(U_mod, J):
    mod_error = np.linalg.norm(U_mod - V_CKM_TARGET) / np.linalg.norm(V_CKM_TARGET)
    j_error = abs(J - J_TARGET) / J_TARGET if J_TARGET != 0 else abs(J)
    return mod_error + 0.1 * j_error

def get_representation_and_u1(M, u1_angle_scale=U1_ANGLE_SCALE):
    try:
        # Get SL(2,C) matrices via polished_holonomy
        rho = M.polished_holonomy()
        gen_names = rho.generators()
        if len(gen_names) < 2:
            print(f"  Too few generators: {len(gen_names)}")
            return None, None
        m1 = np.array(rho.SL2C(gen_names[0]), dtype=complex)
        m2 = np.array(rho.SL2C(gen_names[1]), dtype=complex)
        matrices = [m1, m2]
        if len(gen_names) < 3:
            m3 = m1 @ m2
            matrices.append(m3)
        else:
            m3 = np.array(rho.SL2C(gen_names[2]), dtype=complex)
            matrices.append(m3)

        # Get fundamental group and abelianization
        G = M.fundamental_group()
        abel = G.abelianization()
        if abel:
            free_rank = len(abel[0])
        else:
            free_rank = 0

        # Deterministic random phases based on manifold volume
        random.seed(int(M.volume()*1000) if M.volume() else 0)
        free_angles = [random.uniform(-u1_angle_scale, u1_angle_scale) for _ in range(free_rank)]

        # Convert abelianization vectors to numpy arrays
        img = [np.array(abel[i]) for i in range(len(abel))]

        # Vectors for the three elements we use
        img_elem = []
        img_elem.append(img[0])
        img_elem.append(img[1])
        if len(gen_names) < 3:
            img_elem.append(img[0] + img[1])
        else:
            if len(img) > 2:
                img_elem.append(img[2])
            else:
                img_elem.append(img[0] + img[1])

        # Compute U(1) phases
        u1_phases = []
        for vec in img_elem:
            angle = np.dot(vec, free_angles)
            u1_phases.append(np.exp(1j * angle))

        return matrices, u1_phases

    except Exception as e:
        print(f"  Error in get_representation_and_u1: {e}")
        return None, None

def scan_manifolds(limit=50, sigma=0.5, outfile='scan_results.csv', u1_scale=U1_ANGLE_SCALE):
    results = []
    with open(outfile, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'name', 'volume', 'fitness', 'jarlskog',
                         'u11', 'u12', 'u13', 'u21', 'u22', 'u23', 'u31', 'u32', 'u33'])

    for i in range(limit):
        M = snappy.OrientableClosedCensus[i]
        name = M.name()
        try:
            vol = float(M.volume())
        except:
            vol = -1.0
        print(f"Processing {name} (vol={vol:.4f})...")

        try:
            matrices, u1_phases = get_representation_and_u1(M, u1_scale)
            if matrices is None or u1_phases is None or len(matrices) < 3:
                print("  Could not get representation or U(1) phases, skipping.")
                continue

            points = []
            for g in matrices:
                try:
                    pt = fixed_point(g)
                    points.append(pt)
                except:
                    points.append(np.inf)

            U = build_mixing_matrix(points, u1_phases, sigma)
            if U is None:
                print("  Failed to build mixing matrix.")
                continue

            U_mod = np.abs(U)
            J = jarlskog(U)
            score = fitness(U_mod, J)

            results.append({
                'name': name,
                'volume': vol,
                'score': score,
                'J': J,
                'U_mod': U_mod
            })

            with open(outfile, 'a', newline='') as f:
                writer = csv.writer(f)
                flat = [i, name, f"{vol:.6f}", f"{score:.6f}", f"{J:.2e}"]
                flat.extend([f"{x:.4f}" for row in U_mod for x in row])
                writer.writerow(flat)

            print(f"  score={score:.4f}, J={J:.2e}")
        except Exception as e:
            print(f"  Error: {e}")
            continue

    results.sort(key=lambda x: x['score'])
    return results

if __name__ == "__main__":
    sigma = 0.5
    outfile = 'scan_results.csv'
    u1_scale = U1_ANGLE_SCALE
    if len(sys.argv) > 1:
        try:
            sigma = float(sys.argv[1])
        except:
            print(f"Invalid sigma value '{sys.argv[1]}', using default {sigma}")
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    if len(sys.argv) > 3:
        try:
            u1_scale = float(sys.argv[3])
        except:
            print(f"Invalid U(1) scale '{sys.argv[3]}', using default {u1_scale}")
    print(f"Scanning first 50 closed manifolds with sigma={sigma}, U(1) scale={u1_scale}...\n")
    results = scan_manifolds(limit=50, sigma=sigma, outfile=outfile, u1_scale=u1_scale)
    print("\nTop 5 candidates:")
    for r in results[:5]:
        print(f"{r['name']:20s} vol {r['volume']:.4f}  score {r['score']:.4f}  J {r['J']:.2e}")
