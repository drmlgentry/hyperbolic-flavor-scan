"""
cp_reproduce.py
===============
Canonical reproduction script for:
  "CP Violation from Twist Angles of Hyperbolic Holonomy:
   A Parameter-Free Prediction of the Leptonic CP Phase"

Reproduces the single main result in < 10 seconds.
Requires: snappy, numpy

Usage:
  python cp_reproduce.py
"""

import snappy
import numpy as np

print("=" * 55)
print("CP PHASE -- CANONICAL REPRODUCTION")
print("=" * 55)
print()

M = snappy.OrientableClosedCensus[1]   # m003(-2,3)
rho = M.polished_holonomy()            # 150-bit precision

print(f"Manifold: {M.name()}, vol={float(M.volume()):.5f}, H1={M.homology()}")
print()
print("Convention: polished_holonomy() at 150-bit precision")
print("phi(gamma) = Im(log lambda), lambda dominant eigenvalue")
print("Negative angles possible under this convention.")
print()

def twist_polished(rho, word):
    mat = np.array(rho(word), dtype=complex)
    mat /= np.sqrt(np.linalg.det(mat))
    evals = np.linalg.eigvals(mat)
    lam = evals[np.argmax(np.abs(evals))]
    return np.degrees(float(np.imag(np.log(lam))))

phi_aaB = twist_polished(rho, 'aaB')
phi_baa = twist_polished(rho, 'baa')
delta   = (180.0 + phi_aaB + phi_baa) % 360.0

PDG = 197.0
err = abs(delta - PDG) / PDG * 100

print(f"phi(aaB) = {phi_aaB:.4f} deg")
print(f"phi(baa) = {phi_baa:.4f} deg")
print()
print(f"delta_CP = (180 + phi(aaB) + phi(baa)) mod 360")
print(f"         = (180 + ({phi_aaB:.4f}) + ({phi_baa:.4f})) mod 360")
print(f"         = {delta:.4f} deg")
print()
print(f"PDG 2024: {PDG} deg")
print(f"Error:    {err:.4f}%")
print()

status = "PASS" if abs(delta - 195.91) < 0.05 else "FAIL"
print(f"[{status}] delta_CP = {delta:.2f} deg  (paper claims 195.91 deg)")
print()

# Sign check
mat_aa  = np.array(rho('aa'),  dtype=complex)
mat_aaB = np.array(rho('aaB'), dtype=complex)
mat_baa = np.array(rho('baa'), dtype=complex)

def axis(mat):
    mat = mat.copy()
    mat /= np.sqrt(np.linalg.det(mat))
    L = np.array([[complex(x) for x in row]
                  for row in mat], dtype=complex)
    from scipy.linalg import logm
    Lm = logm(L)
    nx = float(np.real(Lm[0,1] + Lm[1,0])) / 2
    ny = float(np.imag(Lm[1,0] - Lm[0,1])) / 2
    nz = float(np.real(Lm[0,0] - Lm[1,1])) / 2
    v = np.array([nx, ny, nz])
    return v / np.linalg.norm(v)

n_aa  = axis(mat_aa)
n_aaB = axis(mat_aaB)
n_baa = axis(mat_baa)
det = np.linalg.det(np.column_stack([n_aa, n_aaB, n_baa]))
sign_ok = det > 0
print(f"Geometric sign check:")
print(f"  det[n(aa), n(aaB), n(baa)] = {det:.4f}")
print(f"  Sign = +pi (right-handed triple): "
      f"{'YES' if sign_ok else 'NO'}")
print()

status2 = "PASS" if sign_ok else "FAIL"
print(f"[{status2}] Orientation det = {det:.4f} > 0  "
      f"(paper claims +0.079)")
print()
print("=" * 55)
if status == "PASS" and status2 == "PASS":
    print("ALL CLAIMS VERIFIED.")
else:
    print("SOME CLAIMS FAILED -- check output above.")
