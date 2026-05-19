import snappy
import numpy as np
from itertools import product as iprod

print("=== HFG-CELESTIAL HOLOGRAPHY CONNECTION COMPUTATIONS ===")
print()

# ── MANIFOLDS ─────────────────────────────────────────────────────────────────
M_pmns = snappy.OrientableClosedCensus[1]   # m003(-2,3), PMNS
M_ckm  = snappy.OrientableClosedCensus[43]  # m006(-5,2), CKM

print(f"PMNS manifold: {M_pmns.name()}  vol={float(M_pmns.volume()):.6f}")
print(f"CKM  manifold: {M_ckm.name()}   vol={float(M_ckm.volume()):.6f}")
print()

# ── PART 1: HOLONOMY EIGENVALUES → PRINCIPAL SERIES DELTA ────────────────────
print("=== PART 1: PRINCIPAL SERIES DELTA FROM HOLONOMY EIGENVALUES ===")
print("For each word: eigenvalue lambda -> Delta = 1 + i*lambda (principal series)")
print()

def get_loxodromic_data(rho, word):
    """Returns (ell, phi_rad, lambda_val, Delta_principal, z_fixed) for a word."""
    try:
        mat = rho(word)
        tr = complex(mat[0][0]) + complex(mat[1][1])
        if abs(tr) <= 2.0:
            return None  # not loxodromic
        disc = tr**2 - 4
        lam  = (tr + np.sqrt(disc)) / 2
        # Geodesic length and twist
        log_lam = np.log(lam)
        ell   = 2 * abs(log_lam.real)
        phi_r = 2 * log_lam.imag
        # Principal series parameter: lambda_val = Im(log lam) so Delta = 1+i*lambda_val
        lambda_val = log_lam.imag   # this is triangulation-independent? TBD
        Delta = complex(1, lambda_val)
        # Fixed point on CP^1 (axis): solve (a*z+b)/(c*z+d) = z
        a = complex(mat[0][0]); b = complex(mat[0][1])
        c = complex(mat[1][0]); d = complex(mat[1][1])
        if abs(c) > 1e-10:
            # z = ((a-d) +/- sqrt((a-d)^2 + 4bc)) / (2c)
            disc2 = (a-d)**2 + 4*b*c
            z1 = ((a-d) + np.sqrt(disc2)) / (2*c)
            z2 = ((a-d) - np.sqrt(disc2)) / (2*c)
            # Take the attracting fixed point (|lam| > 1 corresponds to z1)
            z_fixed = z1 if abs(lam) > 1 else z2
        else:
            z_fixed = complex(np.inf)
        return ell, phi_r, lambda_val, Delta, z_fixed
    except Exception as e:
        return None

# PMNS best triple: aa, aaB, baa
pmns_words = ['aa', 'aaB', 'baa']
ckm_words  = ['aaB', 'AbA', 'AAb']  # best CKM triple

print("--- PMNS m003(-2,3) ---")
rho_pmns = M_pmns.polished_holonomy()
pmns_data = {}
for word in pmns_words:
    result = get_loxodromic_data(rho_pmns, word)
    if result:
        ell, phi_r, lv, Delta, z = result
        pmns_data[word] = result
        print(f"  {word:<6}: ell={ell:.6f}  twist={np.degrees(phi_r):.4f}deg"
              f"  lambda={lv:.6f}  Delta=1+{lv:.4f}i  z={z:.4f}")
    else:
        print(f"  {word}: not loxodromic or error")

print()
print("--- CKM m006(-5,2) ---")
rho_ckm = M_ckm.polished_holonomy()
ckm_data = {}
for word in ckm_words:
    result = get_loxodromic_data(rho_ckm, word)
    if result:
        ell, phi_r, lv, Delta, z = result
        ckm_data[word] = result
        print(f"  {word:<6}: ell={ell:.6f}  twist={np.degrees(phi_r):.4f}deg"
              f"  lambda={lv:.6f}  Delta=1+{lv:.4f}i  z={z:.4f}")

# ── PART 2: CELESTIAL 3-POINT FUNCTION AT HOLONOMY AXES ─────────────────────
print()
print("=== PART 2: CELESTIAL 3-POINT FUNCTION AT HOLONOMY AXIS POSITIONS ===")
print("K(z1,z2,z3; D1,D2,D3) = |z12|^{D3-D1-D2} |z23|^{D1-D2-D3} |z31|^{D2-D3-D1}")
print("(kinematic kernel, dynamical coefficient C stripped off)")
print()

def celestial_3pt_kernel(z1, z2, z3, D1, D2, D3):
    """Celestial 3-point kinematic kernel K_{R1,R2,R3}."""
    z12 = abs(z1 - z2); z23 = abs(z2 - z3); z31 = abs(z3 - z1)
    if z12 < 1e-10 or z23 < 1e-10 or z31 < 1e-10:
        return None
    exp12 = D3 - D1 - D2
    exp23 = D1 - D2 - D3
    exp31 = D2 - D3 - D1
    return (z12**exp12) * (z23**exp23) * (z31**exp31)

if len(pmns_data) == 3:
    words_p = list(pmns_data.keys())
    zs = [pmns_data[w][4] for w in words_p]
    Ds = [pmns_data[w][3] for w in words_p]
    z1,z2,z3 = zs; D1,D2,D3 = Ds

    K = celestial_3pt_kernel(z1, z2, z3, D1, D2, D3)
    print(f"PMNS axis positions: z1={z1:.4f}  z2={z2:.4f}  z3={z3:.4f}")
    print(f"PMNS Delta values:   D1={D1}  D2={D2}  D3={D3}")
    if K is not None:
        print(f"Celestial 3pt kernel K = {K:.6f}  |K|={abs(K):.6f}")
    else:
        print("Kernel undefined (axes coincide)")

    # Pairwise separations on S^2
    print()
    print("Pairwise celestial separations |z_ij|:")
    for i,j in [(0,1),(1,2),(0,2)]:
        sep = abs(zs[i]-zs[j])
        print(f"  |z_{words_p[i]}-z_{words_p[j]}| = {sep:.6f}")

    # Compare to HFG Borel parameters lambda_ij ~ sqrt(2)*d_ij
    print()
    print("HFG Borel comparison (lambda_ij = sqrt(2)*d_ij):")
    for i,j in [(0,1),(1,2),(0,2)]:
        d_ij = abs(zs[i]-zs[j])  # celestial separation as proxy for d_ij
        lam_pred = np.sqrt(2) * d_ij
        print(f"  lambda_{words_p[i]}{words_p[j]} predicted = {lam_pred:.6f}")

# ── PART 3: SELBERG TRACE FORMULA PREPARATION ────────────────────────────────
print()
print("=== PART 3: GEODESIC LENGTH SPECTRUM (Selberg trace formula input) ===")
print()

PHI = (1+5**0.5)/2

for label, M in [("PMNS m003(-2,3)", M_pmns), ("CKM m006(-5,2)", M_ckm)]:
    print(f"--- {label} ---")
    try:
        ls = M.length_spectrum(cutoff=3.0)
        print(f"  {'Length':>12}  {'Mult':>6}  {'ell/log(phi)':>14}  {'Lucas?':>8}")
        for g in ls[:15]:
            ell = float(g.length.real)
            ratio = ell / np.log(PHI)
            # Check if ratio is close to an integer (Lucas index)
            nearest_int = round(ratio)
            is_lucas_like = abs(ratio - nearest_int) < 0.05
            lucas_flag = f"~L_{nearest_int}" if is_lucas_like else ""
            print(f"  {ell:>12.6f}  {g.multiplicity:>6}  {ratio:>14.6f}  {lucas_flag:>8}")
    except Exception as e:
        print(f"  Error: {e}")
    print()

# ── PART 4: PLANCHEREL MEASURE AT LUCAS SPECTRUM ─────────────────────────────
print("=== PART 4: PLANCHEREL MEASURE AT LUCAS PRINCIPAL SERIES VALUES ===")
print("rho_J(Delta) = -1/(8*pi^3) * (Delta-J-1)(Delta+J-1)")
print("Evaluated at Delta = 1 + i*lambda_k where lambda_k = log(L_k)/log(phi)")
print()

lucas_primes = [2, 3, 7, 11, 29]
def plancherel(Delta, J=0):
    return -1/(8*np.pi**3) * (Delta-J-1)*(Delta+J-1)

print(f"{'L_k':>6}  {'lambda_k':>10}  {'Delta':>20}  {'|rho_0(Delta)|':>16}")
print("-"*58)
for L in lucas_primes:
    lam_k = np.log(L) / np.log(PHI)
    Delta = complex(1, lam_k)
    rho = plancherel(Delta, J=0)
    print(f"{L:>6}  {lam_k:>10.6f}  {Delta}  {abs(rho):>16.6f}")

print()
print("For reference, Plancherel zeros at Delta=1 (J=0): rho=0")
print("Lucas primes give non-zero Plancherel weights -> legitimate representations")

print()
print("=== DONE ===")
