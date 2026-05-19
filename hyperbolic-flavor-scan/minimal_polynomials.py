"""
minimal_polynomials.py
======================
Compute the minimal polynomials of lambda(a) and lambda(b) over
the trace fields of m003 and m006.

m003: trace field Q(sqrt(-3)), ring Z[omega]
m006: trace field Q(sqrt(17)), ring Z[sqrt(17)]

Goal: understand what algebraic numbers the generator lengths are,
and whether 2*ell_b - ell_a = ln(mb/mc) has a number-theoretic explanation.

Run: conda run -n sage python minimal_polynomials.py
"""
import snappy
import numpy as np
from numpy.polynomial import polynomial as P

M_P = snappy.OrientableClosedCensus[1]   # m003(-2,3)
M_C = snappy.OrientableClosedCensus[43]  # m006(-5,2)

rho_P = M_P.polished_holonomy()
rho_C = M_C.polished_holonomy()

def get_eigenvalue(rho, word):
    """Dominant eigenvalue of holonomy matrix."""
    mat = np.array(rho(word), dtype=complex)
    mat /= np.sqrt(np.linalg.det(mat))
    eigs = np.linalg.eigvals(mat)
    return eigs[np.argmax(np.abs(eigs))]

def get_matrix(rho, word):
    """Normalized holonomy matrix."""
    mat = np.array(rho(word), dtype=complex)
    mat /= np.sqrt(np.linalg.det(mat))
    return mat

def trace(rho, word):
    mat = get_matrix(rho, word)
    return mat[0,0] + mat[1,1]

print("="*65)
print("GENERATOR EIGENVALUES AND TRACES")
print("="*65)
print()

for label, M, rho, tf in [
    ("m003(-2,3)", M_P, rho_P, "Q(sqrt(-3))"),
    ("m006(-5,2)", M_C, rho_C, "Q(sqrt(17))"),
]:
    print(f"── {label}  trace field: {tf} ──────────────────")
    print()
    G = M.fundamental_group()
    gens = list(G.generators())
    print(f"Generators: {gens}")
    print()

    for w in ['a', 'b', 'aa', 'bb', 'ab', 'aaB', 'baa']:
        lam = get_eigenvalue(rho, w)
        tr  = trace(rho, w)
        ell = 2*abs(float(np.real(np.log(lam))))
        phi = np.degrees(float(np.imag(np.log(lam))))
        print(f"  word={w:>5}: lambda={lam:.6f}"
              f"  tr={complex(tr):.6f}"
              f"  ell={ell:.6f}  phi={phi:.4f}deg")
    print()

print("="*65)
print("TRACE FIELD VERIFICATION")
print("="*65)
print()

# m003: trace field Q(sqrt(-3))
# tr(rho(a)) should lie in Q(sqrt(-3))
# tr(rho(a)) = Re + Im*i
# If in Q(sqrt(-3)) = Q(i*sqrt(3)):
# a + b*i*sqrt(3) for rational a,b
print("m003 traces vs Q(sqrt(-3)):")
sqrt3 = np.sqrt(3)
for w in ['a','b','aa','ab','aB','aaB','baa']:
    tr = complex(trace(rho_P, w))
    re, im = float(tr.real), float(tr.imag)
    # If tr = a + b*omega where omega = e^(2pi*i/3) = -1/2 + i*sqrt(3)/2
    # then tr = a + b*(-1/2 + i*sqrt(3)/2) = (a - b/2) + i*(b*sqrt(3)/2)
    # So Im(tr) = b*sqrt(3)/2 => b = 2*Im(tr)/sqrt(3)
    b_coeff = 2*im/sqrt3
    a_coeff = re + b_coeff/2
    print(f"  tr({w:>4}) = {re:>8.5f} + {im:>8.5f}i"
          f"  = {a_coeff:.5f} + {b_coeff:.5f}*omega"
          f"  [a={a_coeff:.4f}, b={b_coeff:.4f}]")

print()
print("m006 traces vs Q(sqrt(17)):")
sqrt17 = np.sqrt(17)
for w in ['a','b','aa','ab','aB','aaB','baa']:
    tr = complex(trace(rho_C, w))
    re, im = float(tr.real), float(tr.imag)
    # If tr in Q(sqrt(17)) then Im(tr) should be near 0
    # and Re(tr) = p + q*sqrt(17) for rational p,q
    # q = coefficient of sqrt(17)
    # Try: Re(tr) = p + q*sqrt(17)
    # At two rational points... need two equations
    # Use: tr(rho(aa)) = 3 - sqrt(17) [known]
    # So p=3, q=-1 for aa
    print(f"  tr({w:>4}) = {re:>8.5f} + {im:>8.5f}i")

print()
# Verify tr(aa) = 3 - sqrt(17)
tr_aa = complex(trace(rho_C, 'aa'))
predicted = 3 - sqrt17
print(f"Verification: tr(aa on m006)")
print(f"  Computed:  {float(tr_aa.real):.6f} + {float(tr_aa.imag):.6f}i")
print(f"  3-sqrt(17) = {predicted:.6f}")
print(f"  Match: {abs(float(tr_aa.real) - predicted) < 0.001}")
print()

print("="*65)
print("EIGENVALUE ALGEBRAIC STRUCTURE")
print("="*65)
print()

# For m003: lambda(a) satisfies lambda^2 - tr(a)*lambda + 1 = 0
# since det=1. So lambda = (tr +/- sqrt(tr^2-4))/2
print("m003: lambda(a) satisfies lambda^2 - tr(a)*lambda + 1 = 0")
tr_a_P = complex(trace(rho_P, 'a'))
disc_a_P = tr_a_P**2 - 4
lam_a_P = get_eigenvalue(rho_P, 'a')
print(f"  tr(a) = {tr_a_P:.6f}")
print(f"  disc  = tr^2-4 = {disc_a_P:.6f}")
print(f"  sqrt(disc) = {np.sqrt(disc_a_P):.6f}")
print(f"  lambda(a) = {lam_a_P:.6f}")
print(f"  (tr + sqrt(disc))/2 = {(tr_a_P + np.sqrt(disc_a_P))/2:.6f}")
print()

print("m003: lambda(b) satisfies lambda^2 - tr(b)*lambda + 1 = 0")
tr_b_P = complex(trace(rho_P, 'b'))
disc_b_P = tr_b_P**2 - 4
lam_b_P = get_eigenvalue(rho_P, 'b')
print(f"  tr(b) = {tr_b_P:.6f}")
print(f"  disc  = tr^2-4 = {disc_b_P:.6f}")
print(f"  sqrt(disc) = {np.sqrt(disc_b_P):.6f}")
print(f"  lambda(b) = {lam_b_P:.6f}")
print()

print("m006: lambda(a) satisfies lambda^2 - tr(a)*lambda + 1 = 0")
tr_a_C = complex(trace(rho_C, 'a'))
disc_a_C = tr_a_C**2 - 4
lam_a_C = get_eigenvalue(rho_C, 'a')
print(f"  tr(a) = {tr_a_C:.6f}")
print(f"  disc  = tr^2-4 = {disc_a_C:.6f}")
print(f"  sqrt(disc) = {np.sqrt(disc_a_C):.6f}")
print(f"  lambda(a) = {lam_a_C:.6f}")
print()

print("m006: lambda(b) satisfies lambda^2 - tr(b)*lambda + 1 = 0")
tr_b_C = complex(trace(rho_C, 'b'))
disc_b_C = tr_b_C**2 - 4
lam_b_C = get_eigenvalue(rho_C, 'b')
print(f"  tr(b) = {tr_b_C:.6f}")
print(f"  disc  = tr^2-4 = {disc_b_C:.6f}")
print(f"  sqrt(disc) = {np.sqrt(disc_b_C):.6f}")
print(f"  lambda(b) = {lam_b_C:.6f}")
print()

print("="*65)
print("COMPLEX LENGTHS AS ALGEBRAIC NUMBERS")
print("="*65)
print()

print("The complex length of gamma is:")
print("  L(gamma) = ell + i*phi = 2*log(lambda)")
print("where lambda is the dominant eigenvalue.")
print()
print("For m003, generator a:")
L_a_P = 2*np.log(lam_a_P)
print(f"  L(a) = 2*log(lambda(a)) = {L_a_P:.6f}")
print(f"  ell(a) = {abs(L_a_P.real):.6f}")
print(f"  phi(a) = {np.degrees(L_a_P.imag):.4f} deg")
print()
print("  lambda(a) = exp(L(a)/2)")
print(f"  |lambda(a)| = exp(ell(a)/2) = {np.exp(abs(L_a_P.real)/2):.6f}")
print()

print("For m003, generator b:")
L_b_P = 2*np.log(lam_b_P)
print(f"  L(b) = 2*log(lambda(b)) = {L_b_P:.6f}")
print(f"  ell(b) = {abs(L_b_P.real):.6f}")
print(f"  phi(b) = {np.degrees(L_b_P.imag):.4f} deg")
print()

print("="*65)
print("THE mb/mc EQUATION")
print("="*65)
print()

ell_a = abs(L_a_P.real)
ell_b = abs(L_b_P.real)
combination = 2*ell_b - ell_a
mb_mc = 4.18/1.27
print(f"Claim: exp(2*ell_b - ell_a) = mb/mc")
print(f"  2*ell_b - ell_a = 2*{ell_b:.6f} - {ell_a:.6f}")
print(f"                  = {combination:.6f}")
print(f"  exp({combination:.6f}) = {np.exp(combination):.6f}")
print(f"  mb/mc (PDG)     = {mb_mc:.6f}")
print(f"  Error           = {abs(np.exp(combination)-mb_mc)/mb_mc*100:.4f}%")
print()
print(f"Equivalently: 2*ell_b - ell_a = ln(mb/mc)")
print(f"  LHS = {combination:.6f}")
print(f"  RHS = {np.log(mb_mc):.6f}")
print(f"  Diff = {abs(combination - np.log(mb_mc)):.6f}")
print()

print("In terms of traces:")
print("  ell_a = 2*Re(log(lambda(a)))")
print("       = 2*Re(log((tr(a) + sqrt(tr(a)^2-4))/2))")
print()
print("This is a transcendental function of tr(a).")
print("tr(a) lies in Q(sqrt(-3)), the trace field.")
print()
print(f"  tr(a) = {tr_a_P:.6f}")
print(f"  tr(b) = {tr_b_P:.6f}")
print()
print("The equation 2*ell_b - ell_a = ln(mb/mc) is equivalent to:")
print("  2*Re(log(lambda_b)) - Re(log(lambda_a)) = ln(mb/mc)/2")
print("  Re(log(lambda_b^2 / lambda_a)) = ln(mb/mc)/2")
print("  |lambda_b^2 / lambda_a| = sqrt(mb/mc)")
lam_ratio = lam_b_P**2 / lam_a_P
print(f"  |lambda_b^2/lambda_a| = {abs(lam_ratio):.6f}")
print(f"  sqrt(mb/mc) = {np.sqrt(mb_mc):.6f}")
print(f"  Error = {abs(abs(lam_ratio)-np.sqrt(mb_mc))/np.sqrt(mb_mc)*100:.4f}%")
print()
print("So the mb/mc claim is: |lambda_b(m003)^2 / lambda_a(m003)| = sqrt(mb/mc)")
print()
print("This is a claim about the absolute values of algebraic numbers")
print("in a number field extension of Q(sqrt(-3)).")
print("Whether this has arithmetic meaning requires understanding")
print("the Galois action on these eigenvalues.")

print()
print("="*65)
print("DEHN SURGERY PERSPECTIVE")
print("="*65)
print()
print("m003(-2,3) is obtained from the cusped manifold m003")
print("by Dehn filling with slope (-2,3).")
print()
print("The generator translation lengths ell_a, ell_b depend on")
print("the filling slope via the hyperbolic Dehn surgery theorem.")
print("For nearby slopes, ell_a and ell_b vary continuously.")
print()
print("Question: what filling slope (p,q) would give")
print("  2*ell_b(p,q) - ell_a(p,q) = ln(mb/mc) exactly?")
print("  Answer: slope (-2,3) apparently does this to 0.01%.")
print()
print("Adjacent slopes for comparison:")
print("(checking if (-2,3) is special or if nearby slopes also work)")

# Check nearby cusped manifold
try:
    m003_cusped = snappy.Manifold("m003")
    for p,q in [(-2,3),(-1,3),(-3,3),(-2,2),(-2,4),(-1,2),(-3,2)]:
        try:
            M_fill = m003_cusped.copy()
            M_fill.dehn_fill((p,q))
            if not M_fill.solution_type() == 'all tetrahedra positively oriented':
                continue
            rho_f = M_fill.polished_holonomy()
            lam_a_f = get_eigenvalue(rho_f, 'a')
            lam_b_f = get_eigenvalue(rho_f, 'b')
            ell_a_f = 2*abs(float(np.real(np.log(lam_a_f))))
            ell_b_f = 2*abs(float(np.real(np.log(lam_b_f))))
            pred = np.exp(2*ell_b_f - ell_a_f)
            err = abs(pred - mb_mc)/mb_mc*100
            flag = " <-- TARGET" if abs(p)==2 and q==3 else ""
            print(f"  slope ({p:>2},{q}): exp(2*ell_b-ell_a)={pred:.4f}"
                  f"  err={err:.3f}%{flag}")
        except:
            print(f"  slope ({p:>2},{q}): non-hyperbolic or error")
except Exception as e:
    print(f"  Could not check cusped manifold: {e}")
