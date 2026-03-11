import snappy
import numpy as np

M = snappy.OrientableClosedCensus[43]  # m006
rho = M.polished_holonomy()

words = ['aaB', 'AbA', 'AAb']

print("Fixed points on CP^1 for optimal triple on m006")
print("=" * 60)

fixed_pts = []

for w in words:
    mat = np.array(rho(w), dtype=complex)
    det = np.linalg.det(mat)
    mat = mat / np.sqrt(det)
    a, b, c, d = mat[0,0], mat[0,1], mat[1,0], mat[1,1]

    # Fixed points: c*z^2 + (d-a)*z - b = 0
    # z = [(a-d) +/- sqrt((a-d)^2 + 4bc)] / 2c
    disc = (a - d)**2 + 4*b*c
    sqrt_disc = np.sqrt(disc)
    z_plus  = ((a - d) + sqrt_disc) / (2*c)
    z_minus = ((a - d) - sqrt_disc) / (2*c)

    # Stereographic projection to S^2: z -> (x,y,z_coord)
    def stereo_to_sphere(z):
        x = 2*z.real / (1 + abs(z)**2)
        y = 2*z.imag / (1 + abs(z)**2)
        zc = (abs(z)**2 - 1) / (1 + abs(z)**2)
        return np.array([x, y, zc])

    sph_plus  = stereo_to_sphere(z_plus)
    sph_minus = stereo_to_sphere(z_minus)

    # Attracting fixed point: eigenvalue with |lam| < 1 corresponds to attractor
    evals, evecs = np.linalg.eig(mat)
    # attractor is eigenvector for eigenvalue with |lam| < 1
    idx_attr = np.argmin(np.abs(evals))
    evec_attr = evecs[:, idx_attr]
    z_attr = evec_attr[0] / evec_attr[1] if abs(evec_attr[1]) > 1e-10 else complex(np.inf)

    fixed_pts.append({'word': w, 'z_attr': z_attr, 'sph': stereo_to_sphere(z_attr)})

    print(f"Word: {w}")
    print(f"  z+ = {z_plus:.6f}   z- = {z_minus:.6f}")
    print(f"  z_attr (eigenvec) = {z_attr:.6f}")
    print(f"  arg(z_attr) = {np.angle(z_attr):.6f} rad  ({np.degrees(np.angle(z_attr)):.3f} deg)")
    print(f"  |z_attr|    = {abs(z_attr):.6f}")
    print(f"  S^2 point   = ({stereo_to_sphere(z_attr)[0]:.4f}, {stereo_to_sphere(z_attr)[1]:.4f}, {stereo_to_sphere(z_attr)[2]:.4f})")
    print()

print("=" * 60)
print("Angle between S^2 fixed points (compare with axis angles):")
for i in range(3):
    for j in range(i+1, 3):
        s1 = fixed_pts[i]['sph']
        s2 = fixed_pts[j]['sph']
        cos_ang = np.clip(np.dot(s1, s2), -1, 1)
        ang = np.degrees(np.arccos(cos_ang))
        print(f"  {words[i]} <-> {words[j]}: {ang:.3f} deg")

print()
print("Phase from arg(z_attr) -- natural CP^1 phase:")
for fp in fixed_pts:
    print(f"  {fp['word']}: arg(z_attr) = {np.degrees(np.angle(fp['z_attr'])):.4f} deg")

print()
print("Phase differences:")
ph = [np.angle(fp['z_attr']) for fp in fixed_pts]
for i in range(3):
    for j in range(i+1, 3):
        dp = ph[i] - ph[j]
        print(f"  phi({words[i]}) - phi({words[j]}) = {np.degrees(dp):.4f} deg")

# Now try inserting these CP^1 phases into the overlap matrix
print()
print("Testing CP^1 phases in overlap matrix...")

from scipy.linalg import qr

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

Delta = 2.08
O = np.zeros((3,3), dtype=complex)
for i in range(3):
    for j in range(3):
        theta = angle_matrix[i,j]
        K = (1 - np.cos(theta))**(-Delta/2) if i != j else 1.0
        phase = np.exp(1j * (ph[i] - ph[j]))
        O[i,j] = K * phase

U, _ = qr(O)

CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

J = np.imag(U[0,1] * U[1,2] * np.conj(U[0,2]) * np.conj(U[1,1]))
fitness = np.linalg.norm(np.abs(U) - CKM, 'fro')

print(f"  |U| =")
print(np.abs(U).round(5))
print(f"  Frobenius fitness: {fitness:.6f}")
print(f"  Jarlskog J = {J:.6e}  (PDG: 3.00e-05)")
