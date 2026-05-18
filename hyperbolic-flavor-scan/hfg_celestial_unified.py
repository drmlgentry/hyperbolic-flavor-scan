import snappy
import numpy as np
from scipy.linalg import logm, qr
from itertools import permutations

PHI = (1 + 5**0.5) / 2

# ── PDG 2024 targets ──────────────────────────────────────────────────────────
CKM_PDG = np.array([
    [0.97373, 0.22443, 0.00382],
    [0.22438, 0.97314, 0.04214],
    [0.00886, 0.04130, 0.99914]
])

theta12 = np.arcsin(np.sqrt(0.307));  s12,c12 = np.sin(theta12),np.cos(theta12)
theta23 = np.arcsin(np.sqrt(0.546));  s23,c23 = np.sin(theta23),np.cos(theta23)
theta13 = np.arcsin(np.sqrt(0.02219));s13,c13 = np.sin(theta13),np.cos(theta13)
delta = np.radians(197.0); eid = np.exp(1j*delta)
PMNS_complex = np.array([
    [c12*c13,                    s12*c13,                   s13*np.exp(-1j*delta)],
    [-s12*c23-c12*s23*s13*eid,   c12*c23-s12*s23*s13*eid,  s23*c13],
    [ s12*s23-c12*c23*s13*eid,  -c12*s23-s12*c23*s13*eid,  c23*c13]
])
PMNS_PDG = np.abs(PMNS_complex)

# ── Manifolds ─────────────────────────────────────────────────────────────────
M_pmns = snappy.OrientableClosedCensus[1]
M_ckm  = snappy.OrientableClosedCensus[43]

# ── Helper: safe length extraction from SnapPy LengthSpectrum entry ───────────
def safe_length(g):
    """SnapPy length spectrum entries may return complex or method objects."""
    try:
        v = g.length
        # If it's callable (a method), call it
        if callable(v):
            v = v()
        return float(complex(v).real)
    except Exception:
        return None

# ── Helper: axis direction from holonomy (logm / Pauli method) ───────────────
def axis_from_holonomy(rho, word):
    mat = rho(word)
    A = np.array([[complex(mat[0][0]), complex(mat[0][1])],
                  [complex(mat[1][0]), complex(mat[1][1])]], dtype=complex)
    A = A / np.sqrt(np.linalg.det(A))
    L = logm(A)
    a,b,c,d = L[0,0],L[0,1],L[1,0],L[1,1]
    x = float(np.real(b+c))/2
    y = float(np.imag(c-b))/2
    z = float(np.real(a-d))/2
    v = np.array([x,y,z])
    nrm = np.linalg.norm(v)
    return v/nrm if nrm > 1e-12 else np.array([1.,0.,0.])

# ── Helper: attracting fixed point on CP^1 ────────────────────────────────────
def attracting_fixed_point(rho, word):
    mat = rho(word)
    a,b = complex(mat[0][0]),complex(mat[0][1])
    c,d = complex(mat[1][0]),complex(mat[1][1])
    if abs(c) < 1e-12:
        return complex(np.inf) if abs(a-d) < 1e-12 else b/(d-a)
    disc = np.sqrt((d-a)**2 + 4*b*c)
    z1 = ((a-d)+disc)/(2*c)
    z2 = ((a-d)-disc)/(2*c)
    tr = a+d
    s  = np.sqrt(tr*tr - 4)
    lam1 = (tr+s)/2; lam2 = (tr-s)/2
    return z1 if abs(lam1) >= abs(lam2) else z2

# ── Helper: principal series Delta from eigenvalue phase ─────────────────────
def principal_delta(rho, word):
    mat = rho(word)
    A = np.array([[complex(mat[0][0]),complex(mat[0][1])],
                  [complex(mat[1][0]),complex(mat[1][1])]], dtype=complex)
    A = A / np.sqrt(np.linalg.det(A))
    tr = np.trace(A)
    s  = np.sqrt(tr*tr - 4)
    lam = (tr+s)/2
    nu = np.imag(np.log(lam))
    return complex(1, nu)

# ── Helper: celestial kinematic 3-point kernel ────────────────────────────────
def celestial_3pt(z1,z2,z3,D1,D2,D3):
    z12,z23,z31 = abs(z1-z2),abs(z2-z3),abs(z3-z1)
    if min(z12,z23,z31) < 1e-12:
        return None
    return (z12**(D3-D1-D2)) * (z23**(D1-D2-D3)) * (z31**(D2-D3-D1))

# ── Helper: Borel PMNS ────────────────────────────────────────────────────────
def pmns_borel(M, words, signs, use_sqrt2=True):
    rho = M.polished_holonomy()
    axes = [axis_from_holonomy(rho,w) for w in words]
    d12 = np.dot(axes[0],axes[1])
    d13 = np.dot(axes[0],axes[2])
    d23 = np.dot(axes[1],axes[2])
    f = np.sqrt(2.) if use_sqrt2 else 1.
    L = np.array([[1., 0., 0.],
                  [signs[0]*f*d12, 1., 0.],
                  [signs[1]*f*d13, signs[2]*f*d23, 1.]])
    Q,_ = qr(L)
    for col in range(3):
        if Q[0,col] < 0: Q[:,col] = -Q[:,col]
    return np.abs(Q)

# ── Helper: Gaussian CKM ─────────────────────────────────────────────────────
def ckm_gaussian(M, words, sigma):
    rho = M.polished_holonomy()
    axes = [axis_from_holonomy(rho,w) for w in words]
    theta = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            cos_t = np.clip(np.dot(axes[i],axes[j]),-1.,1.)
            theta[i,j] = np.arccos(abs(cos_t))
    O = np.exp(-theta**2/(2*sigma**2))
    Q,_ = qr(O)
    for col in range(3):
        if Q[0,col] < 0: Q[:,col] = -Q[:,col]
    return np.abs(Q)

# ── Best-permutation fitness ──────────────────────────────────────────────────
def best_fitness(U, target):
    bf = float('inf'); bp = None
    for perm in permutations([0,1,2]):
        f = np.linalg.norm(U[np.ix_(perm,perm)] - target, 'fro')
        if f < bf: bf,bp = f,perm
    return bf, bp

# =============================================================================
print("="*70)
print("PART 1: HFG MIXING MATRICES")
print("="*70)

pmns_words = ["aa","aaB","baa"]; pmns_signs = (1.,-1.,1.)
U_pmns = pmns_borel(M_pmns, pmns_words, pmns_signs)
bf_pmns, bp_pmns = best_fitness(U_pmns, PMNS_PDG)
print(f"\nPMNS fitness vs PDG 2024: {bf_pmns:.6f}  (target 0.005087)")
print("Matrix:")
print(np.array2string(U_pmns[np.ix_(bp_pmns,bp_pmns)], precision=6, suppress_small=True))

ckm_words = ["aaB","AbA","AAb"]; sigma = 0.49
U_ckm = ckm_gaussian(M_ckm, ckm_words, sigma)
bf_ckm, bp_ckm = best_fitness(U_ckm, CKM_PDG)
print(f"\nCKM fitness vs PDG 2024:  {bf_ckm:.6f}  (target 0.016482)")
print("Matrix:")
print(np.array2string(U_ckm[np.ix_(bp_ckm,bp_ckm)], precision=6, suppress_small=True))

# =============================================================================
print("\n" + "="*70)
print("PART 2: CELESTIAL THREE-POINT KERNEL AT HOLONOMY FIXED POINTS")
print("="*70)

rho_pmns = M_pmns.polished_holonomy()
rho_ckm  = M_ckm.polished_holonomy()

print("\n--- PMNS m003(-2,3) ---")
pmns_z = [attracting_fixed_point(rho_pmns,w) for w in pmns_words]
pmns_D = [principal_delta(rho_pmns,w) for w in pmns_words]
for w,z,D in zip(pmns_words,pmns_z,pmns_D):
    print(f"  {w:<6}: z={z:.5f}  Delta={D:.5f}")

K_pmns = celestial_3pt(*pmns_z, *pmns_D)
print(f"\n  3pt kernel K = {K_pmns:.5f}  |K| = {abs(K_pmns):.5f}")

# Pairwise separations and Borel comparison
print("\n  Pairwise |z_ij| vs sqrt(2)*|z_ij| (Borel lambda_ij prediction):")
pairs = [(0,1,pmns_words[0],pmns_words[1]),
         (1,2,pmns_words[1],pmns_words[2]),
         (0,2,pmns_words[0],pmns_words[2])]
for i,j,wi,wj in pairs:
    d = abs(pmns_z[i]-pmns_z[j])
    print(f"  |z_{wi}-z_{wj}| = {d:.6f}   lambda_ij(pred) = {np.sqrt(2)*d:.6f}")

# Cross-check: what does sqrt(2)*d_ij give for the Borel matrix?
print()
rho_p = M_pmns.polished_holonomy()
axes_p = [axis_from_holonomy(rho_p,w) for w in pmns_words]
d12 = np.dot(axes_p[0],axes_p[1])
d13 = np.dot(axes_p[0],axes_p[2])
d23 = np.dot(axes_p[1],axes_p[2])
print(f"  Borel dot products (S3 axes): d12={d12:.6f} d13={d13:.6f} d23={d23:.6f}")
print(f"  CP1 separations:              d12={abs(pmns_z[0]-pmns_z[1]):.6f} "
      f"d13={abs(pmns_z[0]-pmns_z[2]):.6f} d23={abs(pmns_z[1]-pmns_z[2]):.6f}")
print(f"  Ratio CP1/S3: "
      f"{abs(pmns_z[0]-pmns_z[1])/abs(d12):.4f}  "
      f"{abs(pmns_z[0]-pmns_z[2])/abs(d13):.4f}  "
      f"{abs(pmns_z[1]-pmns_z[2])/abs(d23):.4f}")

print("\n--- CKM m006(-5,2) ---")
ckm_z = [attracting_fixed_point(rho_ckm,w) for w in ckm_words]
ckm_D = [principal_delta(rho_ckm,w) for w in ckm_words]
for w,z,D in zip(ckm_words,ckm_z,ckm_D):
    print(f"  {w:<6}: z={z:.5f}  Delta={D:.5f}")

K_ckm = celestial_3pt(*ckm_z, *ckm_D)
if K_ckm is not None:
    print(f"\n  3pt kernel K = {K_ckm:.5f}  |K| = {abs(K_ckm):.5f}")
else:
    print("\n  3pt kernel undefined (degenerate axes)")

# =============================================================================
print("\n" + "="*70)
print("PART 3: GEODESIC LENGTH SPECTRUM vs LUCAS LADDER")
print("="*70)

for label, M in [("PMNS m003(-2,3)", M_pmns), ("CKM m006(-5,2)", M_ckm)]:
    print(f"\n--- {label} ---")
    print(f"  {'ell':>10}  {'mult':>5}  {'ell/log(phi)':>13}  {'note':>10}")
    try:
        ls = M.length_spectrum(cutoff=4.0)
        for g in ls[:20]:
            ell = safe_length(g)
            if ell is None:
                continue
            mult = g.multiplicity
            ratio = ell / np.log(PHI)
            nearest = round(ratio)
            diff = abs(ratio - nearest)
            note = f"~{nearest}*logφ" if (nearest > 0 and diff < 0.03) else ""
            print(f"  {ell:>10.6f}  {mult:>5}  {ratio:>13.6f}  {note:>10}")
    except Exception as e:
        print(f"  Error: {e}")

# =============================================================================
print("\n" + "="*70)
print("PART 4: AXIS GEOMETRY — GOLDEN RATIO RATIOS")
print("="*70)

print("\nPMNS CP1 fixed point separations:")
d_pairs = {}
for i,j,wi,wj in pairs:
    d = abs(pmns_z[i]-pmns_z[j])
    d_pairs[(wi,wj)] = d
    print(f"  d({wi},{wj}) = {d:.6f}")

ds = sorted(d_pairs.values())
print(f"\nRatios of separations:")
for k in range(len(ds)):
    for l in range(k+1,len(ds)):
        ratio = ds[l]/ds[k]
        print(f"  {ds[l]:.6f}/{ds[k]:.6f} = {ratio:.6f}  "
              f"  (phi={PHI:.6f}, phi^2={PHI**2:.6f})")

print("\n=== DONE ===")
