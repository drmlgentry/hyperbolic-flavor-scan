"""
check_r2_r5.py
Find minimal polynomials of the zero-regulator roots r2, r5
over K and over QQ. Run with: sage check_r2_r5.py
"""
from sage.all import *
import mpmath
mpmath.mp.dps = 80

QQx = PolynomialRing(QQ, 'x')
x = QQx.gen()
K = NumberField(x**4 - x - 1, 'w')
w = K.gen()
Ky = PolynomialRing(K, 'y')
y = Ky.gen()

# p6 over K
p6 = (y**6 + w*y**5 + w**3*y**4 + (-w**3+w)*y**3
      + (w**2-2*w)*y**2 + (-2*w**2+w)*y + w**2)

print("Factoring p6 over K:")
facs = p6.factor()
print(facs)
print()
print("Factor degrees:", [(f.degree(), e) for f,e in facs])
print()

# Numerical values of all 6 roots
w_c = mpmath.mpc('-0.24812606280262193', '1.0339820609759678')

def eval_Kpoly(poly, z):
    result = mpmath.mpc(0)
    for i, coeff in enumerate(poly.list()):
        c_list = coeff.list() if hasattr(coeff, 'list') else [coeff]
        c_num = sum(float(c_list[j]) * w_c**j
                    for j in range(len(c_list)))
        result += c_num * z**i
    return result

roots_num = {
    'r1': mpmath.mpc('0.67902584193216324', '-0.18724121458557928'),
    'r2': mpmath.mpc('0.61037204230287976', '-0.54634199607828135'),
    'r3': mpmath.mpc('0.61037204230287974', '0.54634199607828136'),
    'r4': mpmath.mpc('-0.92715190473478514', '-0.84674084639038850'),
    'r5': mpmath.mpc('-0.36224597950025779', '1.3318948975388102'),
    'r6': mpmath.mpc('-0.36224597950025789', '-1.3318948975388102'),
}

print("Evaluating each factor at each root:")
for fac, exp in facs:
    print(f"\n  Factor deg {fac.degree()}: {fac}")
    for rname, rval in roots_num.items():
        v = eval_Kpoly(fac, rval)
        if abs(v) < mpmath.mpf('1e-10'):
            print(f"    ZERO at {rname}: |val| = {abs(v):.2e}")

print()
# Try minimal polynomial of r2 over QQ using algdep at high precision
print("Minimal polynomial of r2 over QQ (algdep, degree 4):")
CC200 = ComplexField(200)
r2_sage = CC200('0.61037204230287975961476963702203888103402567500992',
                '-0.54634199607828134618679413134735129695073011449298')
r5_sage = CC200('-0.36224597950025778637479867860858141565652420976888',
                '1.3318948975388102159310351315719220109776536082303')

for deg in [2, 3, 4, 6]:
    try:
        p = algdep(r2_sage, deg, known_bits=150)
        val = abs(p(r2_sage))
        if val < 1e-30:
            print(f"  deg {deg}: {p}")
            print(f"  |p(r2)| = {float(val):.2e}")
            print(f"  factors over QQ: {QQx(p).factor()}")
            break
    except Exception as e:
        print(f"  deg {deg}: error {e}")

print()
print("Minimal polynomial of r5 over QQ:")
for deg in [2, 3, 4, 6]:
    try:
        p = algdep(r5_sage, deg, known_bits=150)
        val = abs(p(r5_sage))
        if val < 1e-30:
            print(f"  deg {deg}: {p}")
            print(f"  |p(r5)| = {float(val):.2e}")
            print(f"  factors over QQ: {QQx(p).factor()}")
            break
    except Exception as e:
        print(f"  deg {deg}: error {e}")

print()
# Check: does u2 = -w^2-1 appear in the splitting of p6 factors?
u2 = -w**2 - 1
print(f"u2 = {u2}, N(u2) = {u2.norm()}")
print()
# Try factoring p6 over K(sqrt(u2))
Ku2y = PolynomialRing(K, 'z2')
z2 = Ku2y.gen()
try:
    L2 = K.extension(z2**2 - u2, 'b')
    L2y = PolynomialRing(L2, 'yL')
    yL = L2y.gen()
    p6_L2 = (yL**6 + L2(w)*yL**5 + L2(w)**3*yL**4
             + (-L2(w)**3+L2(w))*yL**3
             + (L2(w)**2-2*L2(w))*yL**2
             + (-2*L2(w)**2+L2(w))*yL + L2(w)**2)
    facs_L2 = p6_L2.factor()
    degs = sorted([f.degree() for f,e in facs_L2])
    print(f"p6 over K(sqrt(u2)): degrees {degs}")
    if degs != [6]:
        print("  SPLITS! u2 is relevant.")
    else:
        print("  No splitting. u2 not the key.")
except Exception as e:
    print(f"K(sqrt(u2)) error: {e}")
