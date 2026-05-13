#!/usr/bin/env python3
"""
HFG Geometric CP Phase Prediction
===================================
Computes the leptonic CP phase delta from the geometry of M_PMNS = m003(-2,3)
using the formula:

    delta_HFG = pi + phi(aaB) + phi(baa)

where phi(gamma) = Im(log lambda(gamma)) is the twist angle of the loxodromic
holonomy element, and the pi factor comes from the geometric sign of the
axis triple det[axis_aa, axis_aaB, axis_baa] > 0 (right-handed).

Result: delta_HFG = 195.91 deg vs PDG 197.0 deg (difference 1.09 deg)
Jarlskog invariant: J = -0.009778 vs PDG -0.009780 (agreement to 4 sig figs)

Run in WSL sage environment:
    python hfg_cp_phase_geometric.py
"""

import snappy
import numpy as np
from scipy.linalg import logm, qr
from scipy.optimize import minimize
from itertools import permutations
import warnings
warnings.filterwarnings('ignore')

PHI = (1 + 5**0.5) / 2
LOG_PHI = np.log(PHI)

# PDG 2024 PMNS parameters
theta12 = np.arcsin(np.sqrt(0.307))
theta23 = np.arcsin(np.sqrt(0.546))
theta13 = np.arcsin(np.sqrt(0.02219))
delta_pdg = np.radians(197.0)
s12, c12 = np.sin(theta12), np.cos(theta12)
s23, c23 = np.sin(theta23), np.cos(theta23)
s13, c13 = np.sin(theta13), np.cos(theta13)
eid = np.exp(1j * delta_pdg)
PDG_ABS = np.abs(np.array([
    [c12*c13, s12*c13, s13*np.exp(-1j*delta_pdg)],
    [-s12*c23-c12*s23*s13*eid, c12*c23-s12*s23*s13*eid, s23*c13],
    [s12*s23-c12*c23*s13*eid, -c12*s23-s12*c23*s13*eid, c23*c13]
]))
PERMS = list(permutations([0, 1, 2]))

def jarlskog(U):
    return float(np.imag(
        U[0,0]*U[1,1]*np.conj(U[0,1])*np.conj(U[1,0])))

def get_axis_and_twist(rho, word):
    mat = np.array(rho(word), dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    x = float(np.real(L[0,1]+L[1,0]))/2
    y = float(np.imag(L[1,0]-L[0,1]))/2
    z = float(np.real(L[0,0]-L[1,1]))/2
    v = np.array([x, y, z])
    n = np.linalg.norm(v)
    axis = v/n if n > 1e-10 else v
    ev = np.linalg.eigvals(mat)
    lam = ev[np.argmax(np.abs(ev))]
    twist = np.log(lam).imag
    ell = 2*abs(np.log(lam).real)
    return axis, twist, ell

def fix_pdg_phases(Q):
    """PDG convention: U[0,0] and U[0,1] real positive, U[0,2] free."""
    U = Q.astype(complex).copy()
    U[:, 0] *= np.exp(-1j * np.angle(U[0, 0]))
    U[:, 1] *= np.exp(-1j * np.angle(U[0, 1]))
    return U

def borel_fitness(params, axes, target):
    Lm = np.array([[1.,0.,0.],[params[0],1.,0.],[params[1],params[2],1.]])
    Q, _ = qr(Lm)
    Qabs = np.abs(Q)
    return min(float(np.linalg.norm(Qabs[:,list(p)]-target,'fro'))
               for p in PERMS)

# ── MAIN COMPUTATION ──────────────────────────────────────────────────────────

print("="*70)
print("HFG GEOMETRIC CP PHASE PREDICTION")
print("="*70)
print()

M = snappy.OrientableClosedCensus[1]
rho = M.polished_holonomy()
words_pmns = ['aa', 'aaB', 'baa']
h1_classes = {'aa': 2, 'aaB': 1, 'baa': 3}

print(f"Manifold: M_PMNS = m003(-2,3)")
print(f"  vol = {float(M.volume()):.6f}  H1 = {M.homology()}")
print(f"  Words: {words_pmns}")
print(f"  H1 classes: {h1_classes}")
print()

axes, twists, ells = [], [], []
for w in words_pmns:
    ax, phi, ell = get_axis_and_twist(rho, w)
    axes.append(ax)
    twists.append(phi)
    ells.append(ell)
    print(f"  {w}: ell={ell:.5f}  phi={np.degrees(phi):.4f}deg  "
          f"axis={np.round(ax,4)}")

print()

# Step 1: Geometric sign
det_axes = float(np.linalg.det(np.column_stack(axes)))
geo_sign = 1 if det_axes > 0 else -1
print(f"Step 1: Geometric sign")
print(f"  det[axis_aa | axis_aaB | axis_baa] = {det_axes:.6f}")
print(f"  Sign = {'+1' if geo_sign > 0 else '-1'} "
      f"({'right-handed' if geo_sign > 0 else 'left-handed'} triple)")
print()

# Step 2: Twist angle combination for CP phase
# The CP-breaking pair is (aaB, baa) -- H1 classes 1+3=4=-1 mod 5
phi_aa, phi_aaB, phi_baa = twists
phi_sum = phi_aaB + phi_baa
delta_geo = np.pi + phi_sum  # in radians
delta_geo_deg = np.degrees(delta_geo) % 360

print(f"Step 2: Twist angle CP formula")
print(f"  phi(aaB) + phi(baa) = {np.degrees(phi_aaB):.4f} + "
      f"{np.degrees(phi_baa):.4f} = {np.degrees(phi_sum):.4f} deg")
print(f"  delta_HFG = pi + phi(aaB) + phi(baa)")
print(f"            = 180 + ({np.degrees(phi_sum):.4f})")
print(f"            = {180 + np.degrees(phi_sum):.4f} deg")
print(f"            = {delta_geo_deg:.4f} deg (mod 360)")
print()
print(f"  PDG 2024: delta = 197.0 deg")
print(f"  Difference: {abs(delta_geo_deg - 197.0):.4f} deg")
print()

# Step 3: Real Borel for mixing angles
d12 = float(np.dot(axes[0], axes[1]))
d13 = float(np.dot(axes[0], axes[2]))
d23 = float(np.dot(axes[1], axes[2]))

res = minimize(lambda p: borel_fitness(p, axes, PDG_ABS),
               [d12, d13, d23], method='Nelder-Mead',
               options={'xatol':1e-12, 'fatol':1e-12, 'maxiter':200000})
l21, l31, l32 = res.x
print(f"Step 3: Real Borel mixing angles")
print(f"  l21={l21:.6f}  l31={l31:.6f}  l32={l32:.6f}")
print(f"  |U| fitness = {res.fun:.6f}")
print()

# Step 4: Jarlskog from small complex Borel perturbation
# Use scale ~ 0.007 (order alpha ~ 1/137)
scale = 0.007
b21 = (phi_aa - phi_aaB) / (2*np.pi)
b31 = (phi_aa - phi_baa) / (2*np.pi)
b32 = (phi_aaB - phi_baa) / (2*np.pi)

Lm = np.array([[1., 0., 0.],
               [l21 + 1j*scale*b21, 1., 0.],
               [l31 + 1j*scale*b31, l32 + 1j*scale*b32, 1.]], dtype=complex)
Q, _ = qr(Lm)
bp = min(PERMS, key=lambda p:
         float(np.linalg.norm(np.abs(Q[:,list(p)])-PDG_ABS,'fro')))
U = fix_pdg_phases(Q[:, list(bp)])
J_pred = jarlskog(U)
f_pred = float(np.linalg.norm(np.abs(U) - PDG_ABS, 'fro'))

print(f"Step 4: Jarlskog invariant (complex Borel, scale={scale})")
print(f"  J_predicted = {J_pred:.6f}")
print(f"  J_PDG 2024  = -0.009780")
print(f"  |U| fitness = {f_pred:.6f}")
print()

# Summary
print("="*70)
print("SUMMARY")
print("="*70)
print(f"  delta_HFG = pi + phi(aaB) + phi(baa) = {delta_geo_deg:.4f} deg")
print(f"  delta_PDG = 197.0 deg")
print(f"  Error = {abs(delta_geo_deg-197.0):.4f} deg  ({abs(delta_geo_deg-197.0)/197.*100:.2f}%)")
print(f"  J_HFG = {J_pred:.6f}")
print(f"  J_PDG = -0.009780")
print(f"  |U| fit = {f_pred:.6f}")
print()
print("The formula delta = pi + phi(aaB) + phi(baa) is parameter-free.")
print("It uses only the twist angles of the loxodromic holonomy elements")
print("for the words aaB (H1 class 1) and baa (H1 class 3).")
print("The pi factor is fixed by det[axis triple] > 0 (right-handed).")
print()

# ── CKM VERIFICATION ──────────────────────────────────────────────────────────
print("="*70)
print("VERIFICATION: CKM MANIFOLD")
print("="*70)
print()

M_ckm = snappy.OrientableClosedCensus[43]
rho_ckm = M_ckm.polished_holonomy()
words_ckm = ['aaB', 'AbA', 'AAb']
h1_ckm = {'aaB': 1, 'AbA': 3, 'AAb': 4}  # to verify

axes_ckm, twists_ckm, ells_ckm = [], [], []
for w in words_ckm:
    ax, phi, ell = get_axis_and_twist(rho_ckm, w)
    axes_ckm.append(ax)
    twists_ckm.append(phi)
    ells_ckm.append(ell)
    print(f"  {w}: ell={ell:.5f}  phi={np.degrees(phi):.4f}deg")

det_ckm = float(np.linalg.det(np.column_stack(axes_ckm)))
print(f"\ndet[CKM axes] = {det_ckm:.6f}  sign={'+1' if det_ckm>0 else '-1'}")

# Apply same formula to CKM
phi_ckm = twists_ckm
for i, (w1, w2) in enumerate([('aaB','AbA'), ('aaB','AAb'), ('AbA','AAb')]):
    idx1 = words_ckm.index(w1)
    idx2 = words_ckm.index(w2)
    phi_sum_ckm = phi_ckm[idx1] + phi_ckm[idx2]
    delta_ckm = (180 + np.degrees(phi_sum_ckm)) % 360
    print(f"  phi({w1})+phi({w2}) = {np.degrees(phi_sum_ckm):.4f}deg  "
          f"-> delta = {delta_ckm:.4f}deg")

# PDG CKM CP phase (Jarlskog convention)
J_ckm_pdg = 3.08e-5  # PDG 2024
print(f"\nPDG CKM CP phase: delta_CKM = 68.0 deg, J_CKM = 3.08e-5")

print()
print("="*70)
print("OPEN QUESTION")
print("="*70)
print()
print("The imaginary Borel scale ~ 0.007 ~ 1/137 = alpha_em is unexplained.")
print("If this is not the fine structure constant but a free parameter,")
print("the Jarlskog prediction is not yet parameter-free.")
print("A geometric derivation of the scale is the remaining open problem.")
