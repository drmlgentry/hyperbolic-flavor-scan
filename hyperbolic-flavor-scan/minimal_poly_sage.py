#!/usr/bin/env sage
"""
minimal_poly_sage.py
====================
Use SageMath's algebraic number field machinery to:

1. Find the minimal polynomial of lambda_b^2/lambda_a over Q
2. Determine if this algebraic number lies in Q(sqrt(-3))
3. Find the minimal polynomial of exp(2*ell_b - ell_a) = |lambda_b^2/lambda_a|
4. Check whether the generator lengths ell_a, ell_b satisfy
   algebraic relations over the trace field

Key question: is |lambda_b^2/lambda_a| = sqrt(mb/mc) exactly,
or only approximately?

Run: conda run -n sage sage minimal_poly_sage.py
"""

import snappy
import numpy as np
from sage.all import *

print("="*65)
print("ALGEBRAIC NUMBER THEORY OF HFG GENERATOR EIGENVALUES")
print("="*65)
print()

# Load manifolds
M_P = snappy.OrientableClosedCensus[1]
M_C = snappy.OrientableClosedCensus[43]
rho_P = M_P.polished_holonomy()
rho_C = M_C.polished_holonomy()

def get_dominant_eigenvalue(rho, word, prec=150):
    """Get dominant eigenvalue at high precision."""
    mat = matrix(CC, [[rho(word)[i,j] for j in range(2)] for i in range(2)])
    det_sqrt = sqrt(mat.det())
    mat = mat / det_sqrt
    eigs = mat.eigenvalues()
    return max(eigs, key=lambda x: abs(x))

print("Computing eigenvalues at 150-bit precision...")
CC = ComplexField(150)

# Get generator eigenvalues
lam_a_P = get_dominant_eigenvalue(rho_P, 'a')
lam_b_P = get_dominant_eigenvalue(rho_P, 'b')
lam_a_C = get_dominant_eigenvalue(rho_C, 'a')
lam_b_C = get_dominant_eigenvalue(rho_C, 'b')

print(f"m003: lambda_a = {lam_a_P}")
print(f"m003: lambda_b = {lam_b_P}")
print(f"m006: lambda_a = {lam_a_C}")
print(f"m006: lambda_b = {lam_b_C}")
print()

print("="*65)
print("SECTION 1: TRACES AND THEIR MINIMAL POLYNOMIALS")
print("="*65)
print()

# Traces satisfy: tr^2 - tr*lambda - ... wait
# Actually lambda satisfies: lambda^2 - tr*lambda + 1 = 0
# So the minimal polynomial of lambda over Q(tr) is x^2 - tr*x + 1

# First: what is tr(a) on m003 as an algebraic number?
def get_trace(rho, word):
    mat = matrix(CC, [[rho(word)[i,j] for j in range(2)] for i in range(2)])
    det_sqrt = sqrt(mat.det())
    mat = mat / det_sqrt
    return mat.trace()

tr_a_P = get_trace(rho_P, 'a')
tr_b_P = get_trace(rho_P, 'b')
tr_aa_P = get_trace(rho_P, 'aa')
tr_ab_P = get_trace(rho_P, 'ab')
tr_aa_C = get_trace(rho_C, 'aa')

print(f"m003 tr(a)  = {tr_a_P}")
print(f"m003 tr(b)  = {tr_b_P}")
print(f"m003 tr(aa) = {tr_aa_P}")
print(f"m003 tr(ab) = {tr_ab_P}")
print()
print(f"m006 tr(aa) = {tr_aa_C}")
print(f"m006 3-sqrt(17) = {3 - sqrt(CC(17))}")
print()

# Try to find minimal polynomial of tr(a) on m003
# Use the algebraic number recognition
print("Attempting to identify tr(a) on m003 as algebraic number...")
# tr(a) = 0.21945 + 0.91447i
# |tr(a)|^2 = tr(a)*conj(tr(a))
tr_a_norm_sq = abs(tr_a_P)^2
tr_a_re = real(tr_a_P)
tr_a_im = imag(tr_a_P)
print(f"  |tr(a)|^2 = {tr_a_norm_sq}")
print(f"  Re(tr(a)) = {tr_a_re}")
print(f"  Im(tr(a)) = {tr_a_im}")
print()

# tr(a) satisfies tr(a)*conj(tr(a)) = |tr(a)|^2
# and tr(a) + conj(tr(a)) = 2*Re(tr(a))
# So tr(a) is a root of:
# x^2 - 2*Re(tr(a))*x + |tr(a)|^2 = 0
re_val = RealField(150)(tr_a_re)
norm_val = RealField(150)(tr_a_norm_sq)
print(f"tr(a) satisfies: x^2 - {2*re_val}*x + {norm_val} = 0")
print()

# Try algebraic recognition with algdep
try:
    # Use both real and imaginary parts
    RF = RealField(100)
    re_approx = RF(tr_a_re)
    im_approx = RF(tr_a_im)
    
    print("PSLQ/algdep on Re(tr(a)):")
    p_re = algdep(re_approx, 4)
    print(f"  minimal poly of Re(tr(a)): {p_re}")
    print(f"  roots: {p_re.roots(CC)}")
    
    print()
    print("PSLQ/algdep on Im(tr(a)):")
    p_im = algdep(im_approx, 4)
    print(f"  minimal poly of Im(tr(a)): {p_im}")
    print(f"  roots: {p_im.roots(CC)}")
    
    print()
    print("PSLQ/algdep on |tr(a)|^2:")
    p_norm = algdep(RF(tr_a_norm_sq), 4)
    print(f"  minimal poly of |tr(a)|^2: {p_norm}")
    
except Exception as e:
    print(f"  algdep error: {e}")

print()
print("="*65)
print("SECTION 2: EIGENVALUE RATIO lambda_b^2/lambda_a")
print("="*65)
print()

ratio = lam_b_P^2 / lam_a_P
ratio_abs = abs(ratio)
ratio_arg = arg(ratio)

print(f"lambda_b^2/lambda_a = {ratio}")
print(f"|lambda_b^2/lambda_a| = {ratio_abs}")
print(f"arg(lambda_b^2/lambda_a) = {ratio_arg} rad = {ratio_arg*180/pi} deg")
print()

mb_mc = RR(4.18/1.27)
sqrt_mb_mc = sqrt(mb_mc)
print(f"sqrt(mb/mc) = {sqrt_mb_mc}")
print(f"Error: {abs(ratio_abs - sqrt_mb_mc)/sqrt_mb_mc * 100}%")
print()

# Try to find minimal polynomial of |lambda_b^2/lambda_a|
print("algdep on |lambda_b^2/lambda_a| (degree 2,4,6,8):")
RF100 = RealField(100)
val = RF100(ratio_abs)
for deg in [2, 4, 6, 8]:
    try:
        p = algdep(val, deg)
        # Check quality of approximation
        pval = p(CC(ratio_abs))
        if abs(pval) < 1e-10:
            print(f"  degree {deg}: {p}  [|p(val)|={abs(pval):.2e}] GOOD")
        else:
            print(f"  degree {deg}: {p}  [|p(val)|={abs(pval):.2e}]")
    except Exception as e:
        print(f"  degree {deg}: error {e}")

print()
print("algdep on sqrt(mb/mc) for comparison:")
val_target = RF100(sqrt_mb_mc)
for deg in [2, 4]:
    try:
        p = algdep(val_target, deg)
        print(f"  degree {deg}: {p}")
    except Exception as e:
        print(f"  degree {deg}: error {e}")

print()
print("="*65)
print("SECTION 3: GENERATOR LENGTHS AS ALGEBRAIC NUMBERS?")
print("="*65)
print()

ell_a = RealField(100)(2*abs(log(lam_a_P)))
ell_b = RealField(100)(2*abs(log(lam_b_P)))
print(f"ell_a = {ell_a}")
print(f"ell_b = {ell_b}")
print()

# Geodesic lengths are transcendental in general
# But their COMBINATIONS might be algebraic
combo = 2*ell_b - ell_a
print(f"2*ell_b - ell_a = {combo}")
print(f"exp(2*ell_b - ell_a) = {exp(RR(combo))}")
print()

print("algdep on ell_a:")
for deg in [2, 4, 6]:
    try:
        p = algdep(ell_a, deg)
        pval = abs(p(CC(ell_a)))
        flag = " GOOD" if pval < 1e-10 else ""
        print(f"  degree {deg}: {p}  [|p(val)|={pval:.2e}]{flag}")
    except Exception as e:
        print(f"  degree {deg}: error {e}")

print()
print("algdep on ell_b:")
for deg in [2, 4, 6]:
    try:
        p = algdep(ell_b, deg)
        pval = abs(p(CC(ell_b)))
        flag = " GOOD" if pval < 1e-10 else ""
        print(f"  degree {deg}: {p}  [|p(val)|={pval:.2e}]{flag}")
    except Exception as e:
        print(f"  degree {deg}: error {e}")

print()
print("algdep on 2*ell_b - ell_a (the mb/mc exponent):")
for deg in [2, 4, 6, 8]:
    try:
        p = algdep(RF100(combo), deg)
        pval = abs(p(CC(combo)))
        flag = " GOOD" if pval < 1e-10 else ""
        print(f"  degree {deg}: {p}  [|p(val)|={pval:.2e}]{flag}")
    except Exception as e:
        print(f"  degree {deg}: error {e}")

print()
print("="*65)
print("SECTION 4: RATIO ell_b/ell_a")
print("="*65)
print()
ratio_ell = ell_b/ell_a
print(f"ell_b/ell_a = {ratio_ell}")
print(f"Nearest simple fractions:")
for num in range(1,20):
    for den in range(1,20):
        if abs(ratio_ell - num/den) < 0.002:
            print(f"  {num}/{den} = {num/den:.6f}  err={abs(ratio_ell-num/den)*100:.4f}%")

print()
print("algdep on ell_b/ell_a:")
for deg in [2, 4, 6]:
    try:
        p = algdep(ratio_ell, deg)
        pval = abs(p(CC(ratio_ell)))
        flag = " GOOD" if pval < 1e-10 else ""
        print(f"  degree {deg}: {p}  [|p(val)|={pval:.2e}]{flag}")
    except Exception as e:
        print(f"  degree {deg}: error {e}")

print()
print("="*65)
print("SECTION 5: m006 GENERATOR LENGTH CHECK")
print("="*65)
print()

ell_a_C = RealField(100)(2*abs(log(lam_a_C)))
ell_b_C = RealField(100)(2*abs(log(lam_b_C)))
print(f"m006 ell_a = {ell_a_C}")
print(f"m006 ell_b = {ell_b_C}")
print()

mtau_mmu = RR(1.7769/0.10566)
print(f"exp(3*ell_a on m006) = {exp(3*RR(ell_a_C))}")
print(f"mtau/mmu = {mtau_mmu}")
print(f"Error = {abs(exp(3*RR(ell_a_C))-mtau_mmu)/mtau_mmu*100}%")
print()

print("algdep on 3*ell_a (m006):")
val_3a = 3*ell_a_C
for deg in [2, 4, 6]:
    try:
        p = algdep(val_3a, deg)
        pval = abs(p(CC(val_3a)))
        flag = " GOOD" if pval < 1e-10 else ""
        print(f"  degree {deg}: {p}  [|p(val)|={pval:.2e}]{flag}")
    except Exception as e:
        print(f"  degree {deg}: error {e}")

print()
print("="*65)
print("SECTION 6: DIRECT TRACE FIELD ARITHMETIC")  
print("="*65)
print()

# The invariant trace field of m003 is Q(sqrt(-3))
# The trace tr(a*b^{-1}) or other products might lie in this field
# Let's compute tr(gamma) for all short words and check which
# lie in Q(sqrt(-3))

print("Checking which traces lie in Q(sqrt(-3)):")
print("tr in Q(sqrt(-3)) iff tr^2 + 3*Im(tr)^2/... hmm")
print()
print("Better: tr in Q(sqrt(-3)) = Q(omega) where omega=e^{2pi i/3}")
print("iff tr = a + b*omega for RATIONAL a,b")
print("omega = -1/2 + i*sqrt(3)/2")
print()

omega = CC(-1/2 + I*sqrt(CC(3))/2)
sqrt3 = sqrt(RR(3))

def check_in_Qsqrt3(tr, label, tol=0.001):
    """Check if tr = a + b*omega for rational a,b."""
    re_tr = real(tr)
    im_tr = imag(tr)
    # tr = a + b*omega = a + b*(-1/2 + i*sqrt(3)/2)
    # Re(tr) = a - b/2, Im(tr) = b*sqrt(3)/2
    b = 2*im_tr/sqrt3
    a = re_tr + b/2
    # Check if a,b are rational (i.e., close to simple fractions)
    a_frac = RR(a)
    b_frac = RR(b)
    
    # Try to recognize as rationals
    from sage.rings.rational import Rational
    a_rat = None
    b_rat = None
    for num in range(-20, 21):
        for den in range(1, 20):
            if abs(a_frac - num/den) < tol:
                a_rat = (num, den)
            if abs(b_frac - num/den) < tol:
                b_rat = (num, den)
    
    if a_rat and b_rat:
        print(f"  {label}: tr = {a_rat[0]}/{a_rat[1]} + "
              f"{b_rat[0]}/{b_rat[1]}*omega  IN Q(sqrt(-3))")
    else:
        print(f"  {label}: a={float(a_frac):.5f}, b={float(b_frac):.5f}  "
              f"NOT obviously in Q(sqrt(-3)) with small denominators")
    return a_frac, b_frac

for word in ['a','b','aa','bb','ab','aB','aaB','baa','aaB','AbA']:
    try:
        tr = get_trace(rho_P, word)
        check_in_Qsqrt3(tr, f"tr({word} on m003)")
    except Exception as e:
        print(f"  tr({word}): error {e}")

print()
print("Checking m006 traces vs Q(sqrt(17)):")
print("tr in Q(sqrt(17)) iff tr = a + b*sqrt(17) for rational a,b")
print("(Note: representation is complex, but field is real)")
print()

sqrt17 = sqrt(RR(17))
for word in ['aa','bb','ab','aaB']:
    try:
        tr = get_trace(rho_C, word)
        re_tr = float(real(tr))
        # tr(aa) = 3 - sqrt(17) is real embedded non-standardly
        # Try: the norm N(tr) = tr * Galois_conjugate(tr) should be rational
        norm_tr = abs(tr)^2  # This is tr * conj(tr), not Galois norm
        print(f"  tr({word} on m006) = {float(real(tr)):.5f} + "
              f"{float(imag(tr)):.5f}i  |tr|^2={float(norm_tr):.5f}")
    except Exception as e:
        print(f"  tr({word}): error {e}")
