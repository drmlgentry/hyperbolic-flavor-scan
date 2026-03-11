import snappy
import numpy as np
from scipy.linalg import logm, qr
from scipy.optimize import minimize_scalar, differential_evolution

M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()
words = ['aaB', 'AbA', 'AAb']

CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])
J_pdg = 3.00e-5

def matrix_to_axis(mat):
    mat = np.array(mat, dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    x = float(np.real(L[0,1] + L[1,0])) / 2
    y = float(np.imag(L[1,0] - L[0,1])) / 2
    z = float(np.real(L[0,0] - L[1,1])) / 2
    v = np.array([x, y, z])
    return v / np.linalg.norm(v)

axes = [matrix_to_axis(rho(w)) for w in words]
angle_matrix = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        c = np.clip(np.dot(axes[i], axes[j]), -1, 1)
        angle_matrix[i,j] = np.arccos(abs(c))

Delta = 2.08

# Homology images of the three words: aaB->[2,-1], AbA->[0,0], AAb->[-2,1]
# A flat U(1) connection is a character chi: Z/5 -> U(1)
# chi(generator_a) = exp(2*pi*i*k/5), k in {0,1,2,3,4}
# Homology image of word w = n_a*[1,0] + n_b*[0,1] in Z^2 -> Z/5
# We need the abelianization map. For m006, H1=Z/5.
# The abelianization sends a->1, b->? We need to find the map.
# From homology images [2,-1],[0,0],[-2,1]: these are in Z^2 before abelianization.
# Let's parametrize: phi_i = alpha * h_i where h_i is some integer from homology.
# 
# Strategy: treat the three phases phi_1, phi_2, phi_3 as free parameters
# but constrain phi_2 = 0 (AbA has trivial homology, gets no phase)
# and scan (phi_1, phi_3) freely to find best J with acceptable fitness.

def compute_U_and_J(phi1, phi2, phi3, Delta=2.08):
    phases = [phi1, phi2, phi3]
    O = np.zeros((3,3), dtype=complex)
    for i in range(3):
        for j in range(3):
            theta = angle_matrix[i,j]
            K = (1 - np.cos(theta))**(-Delta/2) if i != j else 1.0
            phase = np.exp(1j * (phases[i] - phases[j]))
            O[i,j] = K * phase
    U, _ = qr(O)
    J = np.imag(U[0,1] * U[1,2] * np.conj(U[0,2]) * np.conj(U[1,1]))
    fitness = np.linalg.norm(np.abs(U) - CKM, 'fro')
    return U, J, fitness

# phi_2 = 0 forced by trivial homology of AbA
# Scan phi_1, phi_3 on a grid first
print("Grid scan: phi_1 in [-pi,pi], phi_3 in [-pi,pi], phi_2=0 fixed")
print("Looking for: fitness < 0.05 AND |J| > 1e-6")
print()

best_J = 0
best_result = None
results = []

N = 60
for phi1 in np.linspace(-np.pi, np.pi, N):
    for phi3 in np.linspace(-np.pi, np.pi, N):
        U, J, fitness = compute_U_and_J(phi1, 0.0, phi3)
        results.append((abs(J), fitness, phi1, phi3, J))
        if fitness < 0.05 and abs(J) > abs(best_J):
            best_J = J
            best_result = (phi1, phi3, U, J, fitness)

results.sort(key=lambda x: -x[0])

print("Top 10 by |J| regardless of fitness:")
for r in results[:10]:
    print(f"  |J|={r[0]:.4e}  fitness={r[1]:.4f}  phi1={np.degrees(r[2]):.1f}  phi3={np.degrees(r[3]):.1f}")

print()
good = [(r) for r in results if r[1] < 0.10]
good.sort(key=lambda x: -x[0])
print(f"Results with fitness < 0.10: {len(good)}")
if good:
    print("Top 10 by |J| with fitness < 0.10:")
    for r in good[:10]:
        print(f"  |J|={r[0]:.4e}  fitness={r[1]:.4f}  phi1={np.degrees(r[2]):.1f}  phi3={np.degrees(r[3]):.1f}")

print()
if best_result:
    phi1, phi3, U, J, fitness = best_result
    print(f"Best result with fitness < 0.05:")
    print(f"  phi1 = {np.degrees(phi1):.3f} deg, phi3 = {np.degrees(phi3):.3f} deg")
    print(f"  J = {J:.6e}  (PDG: 3.00e-05)")
    print(f"  fitness = {fitness:.6f}")
    print(f"  |U| =")
    print(np.abs(U).round(5))
else:
    print("No results found with fitness < 0.05")
    print("Relaxing to fitness < 0.10...")
    if good:
        r = good[0]
        phi1, phi3 = r[2], r[3]
        U, J, fitness = compute_U_and_J(phi1, 0.0, phi3)
        print(f"  phi1 = {np.degrees(phi1):.3f} deg, phi3 = {np.degrees(phi3):.3f} deg")
        print(f"  J = {J:.6e}  (PDG: 3.00e-05)")
        print(f"  fitness = {fitness:.6f}")
        print(f"  |U| =")
        print(np.abs(U).round(5))
