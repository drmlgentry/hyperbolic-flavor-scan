"""
verify_sextic_theorem.py
Verifies the complete Sextic-Octic Decomposition Theorem:
  p8 = Q1*Q2*Q3*Q4  (product of all 4 Galois conjugates of q2)
  p6 = Q2*Q3*Q4     (complementary triple)
  Norm_{K/Q}(p6) = p8^3

Run with: sage verify_sextic_theorem.py
All checks must print PASS.
"""
from sage.all import *

print("Verifying Sextic-Octic Decomposition Theorem")
print("="*60)

QQx = PolynomialRing(QQ, 'x')
x = QQx.gen()
QQy = PolynomialRing(QQ, 'y')
y_QQ = QQy.gen()

K = NumberField(x**4 - x - 1, 'w')
w = K.gen()
Ky = PolynomialRing(K, 'y')
y = Ky.gen()

# The three key polynomials
q2 = y**2 - w*y + (-w**3 + w**2 + 1)
p6 = (y**6 + w*y**5 + w**3*y**4 + (-w**3+w)*y**3
      + (w**2-2*w)*y**2 + (-2*w**2+w)*y + w**2)
p8_QQ = y_QQ**8 + y_QQ**6 - 2*y_QQ**5 + y_QQ**4 - y_QQ**3 + 3*y_QQ**2 - 3*y_QQ + 1

# CHECK 1: p8 = Norm_{K/Q}(q2)
print("\nCHECK 1: p8 = Norm_{K/Q}(q2)")
Qyw = PolynomialRing(QQ, ['y','w'])
yy, ww = Qyw.gens()
q2_yw = yy**2 - ww*yy + (-ww**3 + ww**2 + 1)
f_w = ww**4 - ww - 1
norm_q2 = QQy(q2_yw.resultant(f_w, ww))
assert norm_q2 == p8_QQ, f"FAIL: Norm(q2) = {norm_q2}"
print("  PASS: Norm_{K/Q}(q2) = p8")

# CHECK 2: p6 irreducible over K
print("\nCHECK 2: p6 irreducible over K")
facs = p6.factor()
assert len(list(facs)) == 1 and list(facs)[0][1] == 1, f"FAIL: p6 factors as {facs}"
print("  PASS: p6 is irreducible over K")

# CHECK 3: Over splitting field L, p8 = Q1*Q2*Q3*Q4
print("\nCHECK 3: p8 = Q1*Q2*Q3*Q4 over L (splitting field of x^4-x-1)")
Lsplit = (x**4 - x - 1).splitting_field('a')
embs = K.embeddings(Lsplit)
assert len(embs) == 4, f"FAIL: expected 4 embeddings, got {len(embs)}"
RLy = PolynomialRing(Lsplit, 'yL')
yL = RLy.gen()
prod_all4 = RLy(1)
Qi = []
for phi in embs:
    wL = phi(w)
    Qi_factor = yL**2 - wL*yL + (-wL**3 + wL**2 + 1)
    Qi.append(Qi_factor)
    prod_all4 *= Qi_factor
p8_L = RLy(Ky(p8_QQ))
assert prod_all4 == p8_L, "FAIL: product of 4 conjugates != p8"
print("  PASS: Q1*Q2*Q3*Q4 = p8 over L")

# CHECK 4: p6 = product of exactly 3 conjugates (one specific combo)
print("\nCHECK 4: p6 = Qi*Qj*Qk for some triple {i,j,k}")
from itertools import combinations
p6_L = RLy(sum(Lsplit(embs[0](c))*yL**i
               for i,c in enumerate(p6.list())))
found_combo = None
for combo in combinations(range(4), 3):
    prod3 = RLy(1)
    for idx in combo:
        prod3 *= Qi[idx]
    if prod3 == p6_L:
        found_combo = combo
        break
assert found_combo is not None, "FAIL: p6 is not a product of 3 conjugates"
missing = [i for i in range(4) if i not in found_combo][0]
print(f"  PASS: p6 = Q{found_combo[0]+1}*Q{found_combo[1]+1}*Q{found_combo[2]+1}")
print(f"  Omitted conjugate: Q{missing+1} (geometric embedding sigma_{missing+1})")

# CHECK 5: Norm_{K/Q}(p6) = p8^3
print("\nCHECK 5: Norm_{K/Q}(p6) = p8^3")
p6_yw = (yy**6 + ww*yy**5 + ww**3*yy**4 + (-ww**3+ww)*yy**3
         + (ww**2-2*ww)*yy**2 + (-2*ww**2+ww)*yy + ww**2)
norm_p6 = QQy(p6_yw.resultant(f_w, ww))
assert norm_p6 == p8_QQ**3, "FAIL: Norm(p6) != p8^3"
print("  PASS: Norm_{K/Q}(p6) = p8^3")

# CHECK 6: p8 * (1/Q_{geometric}) = p6 over L (i.e. p8 = q2 * p6 over L)
print("\nCHECK 6: p8 = Q_{geometric} * p6 over L")
Q_geom = Qi[missing]
prod_check = Q_geom * p6_L
assert prod_check == p8_L, "FAIL: Q_geom * p6 != p8"
print(f"  PASS: Q{missing+1} * p6 = p8 over L")

# CHECK 7: discriminant of p8
print("\nCHECK 7: disc(p8) = 7 * 11 * 283^2")
d = p8_QQ.discriminant()
assert d == 7 * 11 * 283**2, f"FAIL: disc = {d}"
print(f"  PASS: disc(p8) = {d} = 7 * 11 * 283^2")

# SUMMARY
print()
print("="*60)
print("ALL 7 CHECKS PASS")
print("="*60)
print()
print("Theorem (Sextic-Octic Decomposition):")
print("  Let K=Q(w), w^4=w+1, L=splitting field of x^4-x-1.")
print("  Let q2 = y^2 - wy + (-w^3+w^2+1) in K[y].")
print("  Let Qi = sigma_i(q2) be the four Galois conjugates over L.")
print()
print("  Then:")
print("  (1) p8 = Q1*Q2*Q3*Q4 = Norm_{K/Q}(q2)         [the octic]")
print(f"  (2) p6 = Q{found_combo[0]+1}*Q{found_combo[1]+1}*Q{found_combo[2]+1}"
      f"                           [the sextic, omitting Q{missing+1}]")
print("  (3) Norm_{K/Q}(p6) = p8^3                       [the cube identity]")
print(f"  (4) Q{missing+1} * p6 = p8  over L               [master relation]")
print()
print("  disc(p8) = 7 * 11 * 283^2")
print("  Gal(p8/Q) = [2^4]S4, order 384")
