import snappy
import numpy as np
from scipy.linalg import logm, qr

print("=" * 65)
print("FULL VERIFICATION: gentry-hyperbolic-flavor-ckm")
print("=" * 65)

M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()

vol = float(M.volume())
print(f"\n[1] MANIFOLD")
print(f"    Name:    {M.name()}")
print(f"    Volume:  {vol:.5f}  (paper: 2.02885)  {'OK' if abs(vol-2.02885)<0.0001 else 'MISMATCH'}")
print(f"    H1:      {M.homology()}  (paper: Z/5)")

def matrix_to_axis(mat):
    mat = np.array(mat, dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    x = float(np.real(L[0,1] + L[1,0])) / 2
    y = float(np.imag(L[1,0] - L[0,1])) / 2
    z = float(np.real(L[0,0] - L[1,1])) / 2
    v = np.array([x, y, z])
    return v / np.linalg.norm(v)

words = ['aaB', 'AbA', 'AAb']
paper_axes = {
    'aaB': np.array([ 0.237, -0.899,  0.368]),
    'AbA': np.array([ 0.376,  0.911,  0.170]),
    'AAb': np.array([ 0.091,  0.184,  0.979]),
}

print(f"\n[2] AXIS VECTORS")
axes = {}
for w in words:
    mat = np.array(rho(w), dtype=complex)
    mat_n = mat / np.sqrt(np.linalg.det(mat))
    tr = abs(float(np.trace(mat_n).real))
    n = matrix_to_axis(mat)
    axes[w] = n
    dot = abs(np.dot(n, paper_axes[w]))
    match = "OK" if dot > 0.999 else "MISMATCH"
    print(f"    {w}: |tr|={tr:.4f}  axis_dot={dot:.6f}  [{match}]")
    print(f"         recomputed: ({n[0]:.3f}, {n[1]:.3f}, {n[2]:.3f})")
    print(f"         paper:      ({paper_axes[w][0]:.3f}, {paper_axes[w][1]:.3f}, {paper_axes[w][2]:.3f})")

paper_angles = {('aaB','AbA'): 0.8405, ('aaB','AAb'): 1.3523, ('AbA','AAb'): 1.1943}
paper_deg    = {('aaB','AbA'): 48.2,   ('aaB','AAb'): 77.5,   ('AbA','AAb'): 68.4}

print(f"\n[3] PAIRWISE ANGLES")
angle_matrix = np.zeros((3,3))
pairs = [('aaB','AbA'), ('aaB','AAb'), ('AbA','AAb')]
for (w1, w2) in pairs:
    i, j = words.index(w1), words.index(w2)
    c = float(np.clip(abs(np.dot(axes[w1], axes[w2])), -1, 1))
    theta = np.arccos(c)
    angle_matrix[i,j] = angle_matrix[j,i] = theta
    diff = abs(theta - paper_angles[(w1,w2)])
    match = "OK" if diff < 0.005 else "MISMATCH"
    print(f"    {w1}<->{w2}: {theta:.4f} rad / {np.degrees(theta):.2f} deg  "
          f"(paper: {paper_angles[(w1,w2)]:.4f} / {paper_deg[(w1,w2)]:.1f} deg)  [{match}]")

# GAUSSIAN kernel (sigma=0.49) -- what the paper actually uses
sigma = 0.49
print(f"\n[4] OVERLAP MATRIX  (Gaussian, sigma={sigma})")
O = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        theta = angle_matrix[i,j]
        O[i,j] = np.exp(-theta**2 / (2*sigma**2)) if i != j else 1.0
for row in O:
    print("    ", "  ".join(f"{v:.5f}" for v in row))

U, R = qr(O)
for i in range(3):
    if U[i,i] < 0:
        U[:,i] *= -1

CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])
paper_U = np.array([
    [0.9744, 0.2246, 0.0110],
    [0.2238, 0.9734, 0.0487],
    [0.0216, 0.0450, 0.9988]
])
fitness = np.linalg.norm(np.abs(U) - CKM, "fro")

print(f"\n[5] MIXING MATRIX  |U|")
labels = [["Vud","Vus","Vub"],["Vcd","Vcs","Vcb"],["Vtd","Vts","Vtb"]]
all_match = True
for i in range(3):
    for j in range(3):
        val = abs(U[i,j])
        diff_p = abs(val - paper_U[i,j])
        match = "OK" if diff_p < 0.002 else "CHECK"
        if diff_p >= 0.002: all_match = False
        print(f"    {labels[i][j]}: {val:.5f}  (paper {paper_U[i,j]:.4f}, PDG {CKM[i,j]:.5f})  [{match}]")
fit_match = "OK" if abs(fitness-0.01729)<0.001 else "MISMATCH"
print(f"    Frobenius: {fitness:.6f}  (paper 0.017290)  [{fit_match}]")

print(f"\n[6] SIGMA SENSITIVITY  (scanning for minimum)")
best_fit, best_sigma = 1e9, None
for s in np.linspace(0.30, 1.0, 71):
    Ot = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            Ot[i,j] = np.exp(-angle_matrix[i,j]**2/(2*s**2)) if i!=j else 1.0
    Ut, _ = qr(Ot)
    for i in range(3):
        if Ut[i,i] < 0: Ut[:,i] *= -1
    f = np.linalg.norm(np.abs(Ut) - CKM, "fro")
    if f < best_fit:
        best_fit, best_sigma = f, s
match = "OK" if abs(best_sigma-0.49)<0.02 else "CHECK"
print(f"    Best sigma: {best_sigma:.4f}  (paper: 0.49)  [{match}]")
print(f"    Best fitness: {best_fit:.6f}  (paper: 0.017290)")

print(f"\n[7] EQUIVALENT TRIPLES")
equiv = [['aaB','AbA','AAb'],['aaB','aBa','AAb'],['bAA','AbA','AAb'],
         ['aaB','AbA','Baa'],['bAA','AbA','Baa']]
for triple in equiv:
    try:
        axs = [matrix_to_axis(rho(w)) for w in triple]
        am = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                am[i,j] = np.arccos(float(np.clip(abs(np.dot(axs[i],axs[j])),-1,1)))
        Ot = np.array([[np.exp(-am[i,j]**2/(2*sigma**2)) if i!=j else 1.0
                        for j in range(3)] for i in range(3)])
        Ut, _ = qr(Ot)
        for i in range(3):
            if Ut[i,i] < 0: Ut[:,i] *= -1
        ft = np.linalg.norm(np.abs(Ut) - CKM, "fro")
        match = "OK" if abs(ft-0.01729)<0.002 else "DIFFERENT"
        print(f"    {triple}: {ft:.6f}  [{match}]")
    except Exception as e:
        print(f"    {triple}: ERROR {e}")

J = np.imag(U[0,1]*U[1,2]*np.conj(U[0,2])*np.conj(U[1,1]))
print(f"\n[8] JARLSKOG:  J = {J:.4e}  (paper ~0, PDG 3.00e-05)")

ang_sum = sum(np.degrees(angle_matrix[words.index(w1), words.index(w2)])
              for w1,w2 in pairs)
excess = ang_sum - 180.0
print(f"\n[9] SPHERICAL TRIANGLE")
print(f"    Angle sum: {ang_sum:.2f} deg  (paper: 193.1)")
print(f"    Excess:    {excess:.2f} deg  (Cabibbo: ~13 deg)")
match = "OK" if abs(excess-13.0)<1.5 else "CHECK"
print(f"    [{match}]")

print(f"\n{'='*65}")
overall = all_match and fit_match=="OK"
print(f"SUMMARY: {'ALL CHECKS PASSED' if overall else 'REVIEW ITEMS ABOVE'}")
print(f"{'='*65}")