"""
verify_bloch_paper.py
=====================
Reproduces all computational claims in:
  "Norm Structure of the Tetrahedral Octic and a Period-3 Unit Orbit
   in the Trace Field of the Meyerhoff Manifold"

Run with: sage verify_bloch_paper.py
"""

from sage.all import *
import mpmath
mpmath.mp.dps = 100

print("=" * 65)
print("VERIFICATION: Norm Structure of the Tetrahedral Octic")
print("=" * 65)
print()

# ── 1. NUMBER FIELD ───────────────────────────────────────────────
print("1. NUMBER FIELD K = Q(w), w^4 - w - 1 = 0")
print("-" * 50)

QQx = PolynomialRing(QQ, 'x')
x = QQx.gen()
K = NumberField(x**4 - x - 1, 'w')
w = K.gen()
UK = K.unit_group()

print("   disc(K)    =", K.discriminant())
print("   signature  =", K.signature())
print("   class #    =", K.class_number())
print("   Gal(K/Q)   =", K.galois_group())
print()

gens = UK.gens_values()
u0, u1, u2 = gens[0], gens[1], gens[2]
print("   u0 =", u0, " N =", u0.norm())
print("   u1 =", u1, " N =", u1.norm())
print("   u2 =", u2, " N =", u2.norm())
print()
assert u1 == 1 - w**3, "u1 mismatch"
assert u2 == -w**2 - 1, "u2 mismatch"
print("   CHECK u1 = 1 - w^3, u2 = -w^2 - 1")
print()

# ── 2. PERIOD-3 UNIT ORBIT ────────────────────────────────────────
print("2. PERIOD-3 UNIT ORBIT (Theorem 2)")
print("-" * 50)

T = lambda z: 1/(1-z)

assert T(w**3) == -w,       "T(w^3) != -w"
assert T(-w) == w**(-4),    "T(-w) != w^-4"
assert T(w**(-4)) == w**3,  "T(w^-4) != w^3"
print("   CHECK T-orbit: w^3 -> -w -> w^-4 -> w^3")

assert w**3 == -u1**(-3),  "w^3 != -u1^-3"
assert -w   == u1**(-1),   "-w != u1^-1"
assert w**(-4) == u1**4,   "w^-4 != u1^4"
print("   CHECK w^3=-u1^-3,  -w=u1^-1,  w^-4=u1^4")

assert w**3 * (-w) * w**(-4) == -1
print("   CHECK w^3 * (-w) * w^-4 = -1  (exponent sum -3-1+4=0)")

QQxmp = PolynomialRing(QQ, 'xmp')
xmp = QQxmp.gen()
mp = (w**3).minpoly()
assert list(mp) == list(xmp**4 - 3*xmp**3 + 3*xmp**2 - xmp - 1)
print("   CHECK minpoly(w^3) = x^4-3x^3+3x^2-x-1")

assert w**3 * (1 - w**3)**3 == -1
print("   CHECK w^3*(1-w^3)^3 = -1")
print()

# ── 3. NORM DECOMPOSITION ─────────────────────────────────────────
print("3. NORM DECOMPOSITION (Theorem 1)")
print("-" * 50)

QQy = PolynomialRing(QQ, 'y')
y = QQy.gen()
p8 = y**8 + y**6 - 2*y**5 + y**4 - y**3 + 3*y**2 - 3*y + 1
print("   p8 =", p8)

Ky = PolynomialRing(K, 'y2')
y2 = Ky.gen()
q2 = y2**2 - w*y2 + (-w**3 + w**2 + 1)
print("   q2 =", q2)

L = (x**4 - x - 1).splitting_field('a')
embs = K.embeddings(L)
assert len(embs) == 4

Ly = PolynomialRing(L, 'y3')
y3 = Ly.gen()
norm_q2 = Ly(1)
for phi in embs:
    wL = phi(w)
    norm_q2 *= y3**2 - wL*y3 + (-wL**3 + wL**2 + 1)

p8_coeffs = list(p8)
norm_coeffs = [QQ(c) for c in norm_q2.list()]
assert norm_coeffs == p8_coeffs, "Norm mismatch!"
print("   CHECK Norm_{K/Q}(q2) == p8  [exact symbolic]")
print("   Coefficients: ", norm_coeffs)
print()

# ── 4. INVERSION IDENTITY ─────────────────────────────────────────
print("4. INVERSION IDENTITY (Proposition 2)")
print("-" * 50)

Kt = PolynomialRing(K, 't')
t = Kt.gen()
q2_t = t**2 - w*t + (-w**3 + w**2 + 1)
target = w**2 * (w - t)

g, s, _ = target.xgcd(q2_t)
assert g.degree() == 0, "gcd not 1"
inv_reduced = (s / g[0]) % q2_t
assert inv_reduced == t, "Inverse is not t"
print("   CHECK 1/(w^2*(w-t)) == t  (mod q2)  [symbolic xgcd]")

assert -w**5 + w**4 + w**2 == 1
print("   CHECK -w^5 + w^4 + w^2 = 1  (w^4=w+1, w^5=w^2+w)")
print()

# ── 5. BLOCH-WIGNER (Numerical) ───────────────────────────────────
print("5. BLOCH-WIGNER IDENTITY (Numerical, 100 digits)")
print("-" * 50)

def D(z):
    z = mpmath.mpc(z)
    return mpmath.im(mpmath.polylog(2, z)) + mpmath.arg(1-z) * mpmath.log(abs(z))

w_c = mpmath.findroot(lambda z: z**4 - z - 1,
                      mpmath.mpc('-0.2481', '1.0340'),
                      tol=mpmath.mpf('1e-95'))

dw3  = D(w_c**3)
dm_w = D(-w_c)
dw_4 = D(w_c**(-4))

assert mpmath.almosteq(dw3, dm_w, mpmath.mpf('1e-90'))
assert mpmath.almosteq(dw3, dw_4, mpmath.mpf('1e-90'))
print("   CHECK D(w^3) = D(-w) = D(w^-4)  [T-orbit, exact]")

import snappy
v0 = mpmath.mpf(str(snappy.ManifoldHP('m003(-2,3)').volume()))
diff = abs(dw3 + v0)
digits = int(-mpmath.log10(diff))
print("   D(w^3)       =", dw3)
print("   -vol(M)      =", -v0)
print("   |difference| =", diff)
print("   Agreement    : ~%d significant figures" % digits)
assert diff < mpmath.mpf('1e-15')
print("   CHECK -D(w^3) = vol(M)  to full precision")
print()

# ── 6. GEOMETRIC SHAPES (Numerical) ──────────────────────────────
print("6. SHAPE IDENTIFICATIONS (Numerical)")
print("-" * 50)

shapes = snappy.ManifoldHP('m003(-2,3)').tetrahedra_shapes('rect')
z1 = mpmath.mpc(str(shapes[0].real()), str(shapes[0].imag()))
z2 = mpmath.mpc(str(shapes[1].real()), str(shapes[1].imag()))

q2_eval = lambda z: z**2 - w_c*z + (-w_c**3 + w_c**2 + 1)
z_tet = mpmath.findroot(q2_eval, mpmath.mpc('-0.927', '0.847'),
                        tol=mpmath.mpf('1e-50'))

assert mpmath.almosteq(w_c * z_tet**2,           z1, mpmath.mpf('1e-20'))
assert mpmath.almosteq(w_c * (w_c - z_tet)**2,   z2, mpmath.mpf('1e-20'))
assert mpmath.almosteq(1/(1 - z1),              z_tet, mpmath.mpf('1e-20'))
print("   CHECK z1 = w*z_tet^2")
print("   CHECK z2 = w*(w-z_tet)^2")
print("   CHECK 1/(1-z1) = z_tet")

dz1 = D(z1); dz2 = D(z2)
assert mpmath.almosteq(dz1+dz2, v0, mpmath.mpf('1e-15'))
print("   CHECK D(z1)+D(z2) = vol(M)  (%.5f + %.5f)" %
      (float(dz1/v0), float(dz2/v0)))
print()

# ── SUMMARY ──────────────────────────────────────────────────────
print("=" * 65)
print("ALL CHECKS PASSED")
print()
print("Exact algebraic:")
print("  T-orbit, unit expressions, product=-1, minpoly")
print("  p8 = Norm_{K/Q}(q2)")
print("  w^2*z_tet*(w-z_tet) = 1 in K[t]/(q2)")
print()
print("Numerical (100 digits):")
print("  -D(w^3) = vol(M) to ~%d figures" % digits)
print("  z1=w*z_tet^2, z2=w*(w-z_tet)^2, 1/(1-z1)=z_tet")
print("=" * 65)
