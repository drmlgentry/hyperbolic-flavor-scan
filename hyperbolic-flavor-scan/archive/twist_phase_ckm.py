import snappy
import numpy as np
from scipy.linalg import logm

M = snappy.OrientableClosedCensus[43]  # m006
rho = M.polished_holonomy()

words = ['aaB', 'AbA', 'AAb']

print("Loxodromic data for optimal triple on m006")
print("=" * 60)

phases = []
lengths = []

for w in words:
    mat = np.array(rho(w), dtype=complex)
    det = np.linalg.det(mat)
    mat = mat / np.sqrt(det)

    # Eigenvalues of SL(2,C) matrix: lambda and 1/lambda
    evals = np.linalg.eigvals(mat)
    # Take the one with |lambda| >= 1
    lam = evals[0] if abs(evals[0]) >= abs(evals[1]) else evals[1]

    # lambda = exp((ell + i*phi)/2)
    log_lam = np.log(lam)          # principal branch
    ell = 2 * log_lam.real         # translation length
    phi = 2 * log_lam.imag         # twist angle

    # Trace check: tr = lambda + 1/lambda = 2*cosh((ell+i*phi)/2)
    tr = np.trace(mat)

    lengths.append(ell)
    phases.append(phi)

    print(f"Word: {w}")
    print(f"  lambda       = {lam:.6f}")
    print(f"  log(lambda)  = {log_lam:.6f}")
    print(f"  ell (transl) = {ell:.6f}  rad  ({np.degrees(ell):.3f} deg)")
    print(f"  phi (twist)  = {phi:.6f}  rad  ({np.degrees(phi):.3f} deg)")
    print(f"  |tr|         = {abs(tr):.6f}")
    print()

print("=" * 60)
print("Phase differences (enter overlap matrix as exp(i*delta_phi)):")
for i in range(3):
    for j in range(i+1, 3):
        dphi = phases[i] - phases[j]
        print(f"  phi({words[i]}) - phi({words[j]}) = {dphi:.6f} rad  ({np.degrees(dphi):.3f} deg)")

print()
print("Now testing: does inserting twist phases generate J != 0?")
print()

# Build complex overlap matrix with twist phases
def mixing_with_phases(axes_list, angle_matrix, phases_list, Delta=2.08):
    O = np.zeros((3,3), dtype=complex)
    for i in range(3):
        for j in range(3):
            theta = angle_matrix[i,j]
            K = (1 - np.cos(theta))**(-Delta/2) if i != j else 1.0
            phase = np.exp(1j * (phases_list[i] - phases_list[j]))
            O[i,j] = K * phase
    from scipy.linalg import qr
    U, _ = qr(O)
    return U

# Recompute axes
def matrix_to_axis(mat):
    mat = np.array(mat, dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    from scipy.linalg import logm
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

U = mixing_with_phases(axes, angle_matrix, phases)

# Jarlskog invariant
# J = Im(V_us V_cb V_ub* V_cs*)  ->  indices [0,1],[1,2],[0,2]*,[1,1]*
J = np.imag(U[0,1] * U[1,2] * np.conj(U[0,2]) * np.conj(U[1,1]))
print(f"Mixing matrix |U| with twist phases:")
print(np.abs(U).round(4))
print()
print(f"Jarlskog invariant J = {J:.6e}")
print(f"PDG value:         J = 3.00e-05")
print()

CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])
fitness = np.linalg.norm(np.abs(U) - CKM, 'fro')
print(f"Frobenius fitness with phases: {fitness:.6f}")
print(f"Frobenius fitness without phases: 0.017290")
