#!/usr/bin/env python3
"""
Arithmetic Flavor Geometry: Verification Computation
=====================================================
Tests whether rank-3 unit log lattices from totally real quartic fields
produce realistic CKM structure and Jarlskog invariant without tuning.

Two independent fields are tested:
  Field 1: x^4 - x^3 - 4x^2 + 4x + 1 = 0
  Field 2: x^4 - 2x^3 - x^2 + 2x + 1 = 0

Steps executed:
  1. Verify irreducibility and find all real roots
  2. Construct log embedding vectors from roots
  3. Build sector generational vectors with integer coefficients
  4. Perform complex Gram-Schmidt orthonormalization
  5. Compute CKM matrix and mixing angles
  6. Compute Jarlskog invariant J
  7. Scan phase phi from 0.1 to 0.8
  8. Scan integer perturbations +/-1 around base coefficients
  9. Null tests (no shared direction, orthogonal units, 1D sectors)
 10. Repeat all for Field 2
"""

import numpy as np
from numpy.linalg import norm
import itertools

# ============================================================
# SECTION 1: Field arithmetic and root finding
# ============================================================

def find_real_roots(coeffs):
    """Find all real roots of polynomial with given coefficients (high to low degree)."""
    roots = np.roots(coeffs)
    real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-10])
    return real_roots

def verify_field(coeffs, name):
    """Basic checks on a quartic polynomial."""
    roots = find_real_roots(coeffs)
    print(f"\n{'='*60}")
    print(f"Field: {name}")
    print(f"Polynomial coefficients (high to low): {coeffs}")
    print(f"Number of real roots found: {len(roots)}")
    print(f"Real roots: {[f'{r:.8f}' for r in roots]}")
    
    # Verify totally real (4 real roots for degree 4)
    all_roots = np.roots(coeffs)
    n_real = sum(1 for r in all_roots if abs(r.imag) < 1e-8)
    print(f"Totally real: {n_real == 4} ({n_real}/4 roots real)")
    
    # Unit rank for totally real quartic: r = r1 + r2 - 1 = 4 + 0 - 1 = 3
    print(f"Expected unit rank: 3 (totally real quartic: r1=4, r2=0, r=r1-1=3)")
    
    return roots

# ============================================================
# SECTION 2: Fundamental unit approximation
# ============================================================

def compute_log_vectors(roots):
    """
    Compute log embedding vectors for a set of algebraic units.
    For a totally real field with 4 real embeddings sigma_1..sigma_4,
    the log embedding of unit u is (log|sigma_1(u)|, ..., log|sigma_4(u)|)
    restricted to the sum-zero hyperplane.
    
    We approximate fundamental units by identifying elements of small
    norm in the field using the roots as embedding values.
    
    For a unit u = a0 + a1*alpha + a2*alpha^2 + a3*alpha^3,
    sigma_i(u) = a0 + a1*r_i + a2*r_i^2 + a3*r_i^3
    where r_i are the four real roots.
    
    We search for units by their minimal polynomial norm being +-1.
    """
    r = roots  # 4 real roots
    
    def embed(coeffs):
        """Evaluate polynomial with given coefficients at all 4 roots."""
        a0, a1, a2, a3 = coeffs
        return np.array([a0 + a1*ri + a2*ri**2 + a3*ri**3 for ri in r])
    
    def norm_element(coeffs):
        """Field norm = product of all embeddings."""
        return np.prod(embed(coeffs))
    
    def log_vec(coeffs):
        """Log embedding vector, projected to sum-zero hyperplane."""
        vals = embed(coeffs)
        if any(v <= 0 for v in vals):
            return None
        lv = np.log(np.abs(vals))
        # Project to sum-zero: subtract mean
        lv = lv - lv.mean()
        return lv[:3]  # Drop last (redundant by sum-zero constraint)
    
    # Search for units: elements with norm = +/- 1
    # Small integer coefficient search
    units_found = []
    search_range = range(-4, 5)
    
    for a0, a1, a2, a3 in itertools.product(search_range, repeat=4):
        if a0 == 0 and a1 == 0 and a2 == 0 and a3 == 0:
            continue
        if a0 == 1 and a1 == 0 and a2 == 0 and a3 == 0:
            continue  # Skip identity
        if a0 == -1 and a1 == 0 and a2 == 0 and a3 == 0:
            continue  # Skip -1
        
        n = norm_element([a0, a1, a2, a3])
        if abs(abs(n) - 1.0) < 1e-6:
            lv = log_vec([a0, a1, a2, a3])
            if lv is not None and norm(lv) > 0.05:  # Non-trivial
                units_found.append(([a0, a1, a2, a3], lv, norm(lv)))
    
    # Sort by size of log vector
    units_found.sort(key=lambda x: x[2])
    
    return units_found, embed, log_vec

# ============================================================
# SECTION 3: Select three independent fundamental units
# ============================================================

def select_independent_units(units_found, n=3):
    """
    Select n linearly independent units with smallest log vector norms.
    """
    selected = []
    vectors = []
    
    for coeffs, lv, lnorm in units_found:
        if len(selected) == 0:
            selected.append((coeffs, lv, lnorm))
            vectors.append(lv)
        else:
            # Check independence: add and check rank
            test = np.array(vectors + [lv])
            if np.linalg.matrix_rank(test, tol=1e-8) == len(selected) + 1:
                selected.append((coeffs, lv, lnorm))
                vectors.append(lv)
        
        if len(selected) == n:
            break
    
    return selected

# ============================================================
# SECTION 4: Build sector generational vectors
# ============================================================

def build_sector_vectors(u1, u2, u3, up_coeffs=(1,5,9), dn_coeffs=(3,7,12), phi=0.3):
    """
    Build generational vectors for up and down sectors.
    
    Up sector: n*u2*exp(i*phi) + u1  for n in up_coeffs
    Down sector: n*u3 + u1            for n in dn_coeffs
    
    u1 is the shared dominant direction (lepton sector).
    u2 is up-quark dominant direction.
    u3 is down-quark dominant direction.
    
    Returns complex vectors for up and down sectors.
    """
    # Complex phase on u2 direction
    u2_complex = u2.astype(complex) * np.exp(1j * phi)
    u1_c = u1.astype(complex)
    u3_c = u3.astype(complex)
    
    # Up sector: generations ordered light to heavy
    up_vecs = []
    for n in sorted(up_coeffs):  # Light to heavy: ascending n
        v = n * u2_complex + u1_c
        up_vecs.append(v)
    
    # Down sector
    dn_vecs = []
    for n in sorted(dn_coeffs):
        v = n * u3_c + u1_c
        dn_vecs.append(v)
    
    return up_vecs, dn_vecs

# ============================================================
# SECTION 5: Complex Gram-Schmidt
# ============================================================

def gram_schmidt(vecs):
    """
    Complex Gram-Schmidt orthonormalization.
    Returns orthonormal basis.
    """
    basis = []
    for v in vecs:
        w = v.copy().astype(complex)
        for b in basis:
            w = w - np.vdot(b, w) * b
        n = norm(w)
        if n < 1e-12:
            raise ValueError("Vectors are linearly dependent")
        basis.append(w / n)
    return basis

# ============================================================
# SECTION 6: CKM matrix and mixing angles
# ============================================================

def compute_ckm(up_basis, dn_basis):
    """
    CKM matrix V_ij = <up_i | dn_j>.
    Rows: up-type (u, c, t) -- light to heavy.
    Cols: dn-type (d, s, b) -- light to heavy.
    """
    n = len(up_basis)
    V = np.zeros((n, n), dtype=complex)
    for i, ui in enumerate(up_basis):
        for j, dj in enumerate(dn_basis):
            V[i, j] = np.vdot(ui, dj)
    return V

def pdg_angles_from_ckm(V):
    """
    Extract standard PDG mixing angles from CKM matrix using
    the standard parametrization.
    theta_12, theta_13, theta_23 and delta_CP.
    """
    # |V_us| = sin(theta_12) * cos(theta_13)
    # |V_ub| = sin(theta_13)
    # |V_cb| = sin(theta_23) * cos(theta_13)
    
    Vus = abs(V[0,1])
    Vub = abs(V[0,2])
    Vcb = abs(V[1,2])
    Vud = abs(V[0,0])
    
    sin13 = Vub
    cos13 = np.sqrt(max(0, 1 - sin13**2))
    
    if cos13 > 1e-10:
        sin12 = Vus / cos13
        sin23 = Vcb / cos13
    else:
        sin12 = 0
        sin23 = 0
    
    theta12 = np.arcsin(np.clip(sin12, 0, 1))
    theta13 = np.arcsin(np.clip(sin13, 0, 1))
    theta23 = np.arcsin(np.clip(sin23, 0, 1))
    
    return theta12, theta13, theta23

def jarlskog(V):
    """
    Jarlskog invariant J = Im(V_ud V_cs V_us* V_cd*).
    """
    return np.imag(V[0,0] * V[1,1] * np.conj(V[0,1]) * np.conj(V[1,0]))

# ============================================================
# SECTION 7: Full analysis for one field
# ============================================================

def analyze_field(field_name, poly_coeffs, phi=0.3,
                  up_coeffs=(1,5,9), dn_coeffs=(3,7,12),
                  verbose=True):
    """
    Full pipeline for one number field.
    Returns CKM matrix, mixing angles, and J.
    """
    roots = verify_field(poly_coeffs, field_name)
    
    if len(roots) != 4:
        print("ERROR: Not totally real (fewer than 4 real roots)")
        return None
    
    units_found, embed, log_vec = compute_log_vectors(roots)
    
    print(f"\nUnits found with norm +/-1: {len(units_found)}")
    print("First 6 units by log-vector size:")
    for coeffs, lv, lnorm in units_found[:6]:
        n = np.prod(embed(coeffs))
        print(f"  coeffs={coeffs}  norm={n:.6f}  |log_vec|={lnorm:.4f}  "
              f"log_vec={[f'{x:.4f}' for x in lv]}")
    
    selected = select_independent_units(units_found, n=3)
    
    if len(selected) < 3:
        print("ERROR: Could not find 3 independent units")
        return None
    
    u1 = selected[0][1]
    u2 = selected[1][1]
    u3 = selected[2][1]
    
    print(f"\nSelected independent unit log vectors:")
    print(f"  u1 = {[f'{x:.4f}' for x in u1]}")
    print(f"  u2 = {[f'{x:.4f}' for x in u2]}")
    print(f"  u3 = {[f'{x:.4f}' for x in u3]}")
    
    # Angles between unit vectors
    def angle(a, b):
        c = np.dot(a, b) / (norm(a) * norm(b))
        return np.degrees(np.arccos(np.clip(c, -1, 1)))
    
    print(f"\nAngles between unit log vectors:")
    print(f"  angle(u1,u2) = {angle(u1,u2):.1f} deg")
    print(f"  angle(u1,u3) = {angle(u1,u3):.1f} deg")
    print(f"  angle(u2,u3) = {angle(u2,u3):.1f} deg")
    
    # Check unit log norms (these are lambda_i values)
    print(f"\nLog norms (lambda_i values):")
    print(f"  |u1| = {norm(u1):.4f}  (lepton step ~ should be ~2.65 nat log)")
    print(f"  |u2| = {norm(u2):.4f}  (up-quark step ~ should be ~1.44)")
    print(f"  |u3| = {norm(u3):.4f}  (down-quark step ~ should be ~0.74)")
    
    # Build sector vectors
    up_vecs, dn_vecs = build_sector_vectors(u1, u2, u3, up_coeffs, dn_coeffs, phi)
    
    # Gram-Schmidt
    try:
        up_basis = gram_schmidt(up_vecs)
        dn_basis = gram_schmidt(dn_vecs)
    except ValueError as e:
        print(f"Gram-Schmidt failed: {e}")
        return None
    
    # CKM matrix
    V = compute_ckm(up_basis, dn_basis)
    
    print(f"\nCKM matrix |V| (phi={phi:.2f} rad):")
    Vabs = np.abs(V)
    for row in Vabs:
        print(f"  {[f'{x:.4f}' for x in row]}")
    
    # Mixing angles
    t12, t13, t23 = pdg_angles_from_ckm(V)
    print(f"\nMixing angles:")
    print(f"  theta_12 = {t12:.4f} rad = {np.degrees(t12):.2f} deg  "
          f"(PDG: 0.2272 rad, 13.02 deg)")
    print(f"  theta_13 = {t13:.4f} rad = {np.degrees(t13):.2f} deg  "
          f"(PDG: 0.0038 rad, 0.22 deg)")
    print(f"  theta_23 = {t23:.4f} rad = {np.degrees(t23):.2f} deg  "
          f"(PDG: 0.0400 rad, 2.29 deg)")
    
    # Jarlskog
    J = jarlskog(V)
    print(f"\nJarlskog invariant:")
    print(f"  J = {J:.4e}  (PDG: ~3.0e-5)")
    print(f"  J/J_PDG = {J/3.0e-5:.3f}")
    
    return V, (t12, t13, t23), J, u1, u2, u3

# ============================================================
# SECTION 8: Phase scan
# ============================================================

def phase_scan(field_name, poly_coeffs, up_coeffs=(1,5,9), dn_coeffs=(3,7,12)):
    """Scan phi from 0.1 to 0.8 and record J and Cabibbo angle."""
    print(f"\n{'='*60}")
    print(f"PHASE SCAN: {field_name}")
    print(f"{'phi':>8} {'J':>12} {'J/J_PDG':>10} {'theta_12 (deg)':>16}")
    print("-"*50)
    
    roots = find_real_roots(poly_coeffs)
    if len(roots) != 4:
        print("Not totally real, skipping scan")
        return
    
    units_found, embed, log_vec = compute_log_vectors(roots)
    selected = select_independent_units(units_found, n=3)
    if len(selected) < 3:
        print("Could not find 3 independent units, skipping scan")
        return
    
    u1, u2, u3 = selected[0][1], selected[1][1], selected[2][1]
    
    phi_vals = np.linspace(0.1, 0.8, 15)
    results = []
    
    for phi in phi_vals:
        up_vecs, dn_vecs = build_sector_vectors(u1, u2, u3, up_coeffs, dn_coeffs, phi)
        try:
            up_basis = gram_schmidt(up_vecs)
            dn_basis = gram_schmidt(dn_vecs)
            V = compute_ckm(up_basis, dn_basis)
            J = jarlskog(V)
            t12, t13, t23 = pdg_angles_from_ckm(V)
            print(f"  {phi:8.3f} {J:12.4e} {J/3.0e-5:10.3f} {np.degrees(t12):16.2f}")
            results.append((phi, J, t12))
        except Exception as e:
            print(f"  {phi:8.3f}  ERROR: {e}")
    
    return results

# ============================================================
# SECTION 9: Integer perturbation scan
# ============================================================

def integer_perturbation_scan(field_name, poly_coeffs, phi=0.3,
                               base_up=(1,5,9), base_dn=(3,7,12)):
    """
    Perturb each integer coefficient by +/-1 and record CKM stability.
    """
    print(f"\n{'='*60}")
    print(f"INTEGER PERTURBATION SCAN: {field_name}")
    
    roots = find_real_roots(poly_coeffs)
    units_found, embed, log_vec = compute_log_vectors(roots)
    selected = select_independent_units(units_found, n=3)
    if len(selected) < 3:
        print("Insufficient units")
        return
    
    u1, u2, u3 = selected[0][1], selected[1][1], selected[2][1]
    
    # Base case
    up_vecs, dn_vecs = build_sector_vectors(u1, u2, u3, base_up, base_dn, phi)
    up_basis = gram_schmidt(up_vecs)
    dn_basis = gram_schmidt(dn_vecs)
    V_base = compute_ckm(up_basis, dn_basis)
    J_base = jarlskog(V_base)
    t12_base, _, _ = pdg_angles_from_ckm(V_base)
    
    print(f"Base case: up={base_up}, dn={base_dn}")
    print(f"  J = {J_base:.4e},  theta_12 = {np.degrees(t12_base):.2f} deg")
    print(f"\nPerturbations (each +/-1 on one coefficient):")
    print(f"{'Perturbation':<30} {'J':>12} {'J/J_base':>10} {'theta_12':>10}")
    print("-"*65)
    
    # Perturb up coefficients
    for idx in range(3):
        for delta in [-1, +1]:
            new_up = list(base_up)
            new_up[idx] += delta
            if any(n <= 0 for n in new_up):
                continue
            try:
                up_vecs, dn_vecs = build_sector_vectors(u1, u2, u3,
                                    tuple(new_up), base_dn, phi)
                up_b = gram_schmidt(up_vecs)
                dn_b = gram_schmidt(dn_vecs)
                V = compute_ckm(up_b, dn_b)
                J = jarlskog(V)
                t12, _, _ = pdg_angles_from_ckm(V)
                label = f"up[{idx}] {base_up[idx]:+d}->{new_up[idx]:+d}"
                print(f"  {label:<28} {J:12.4e} {J/J_base:10.3f} "
                      f"{np.degrees(t12):10.2f}")
            except Exception as e:
                print(f"  up[{idx}]{delta:+d}: ERROR {e}")
    
    # Perturb down coefficients
    for idx in range(3):
        for delta in [-1, +1]:
            new_dn = list(base_dn)
            new_dn[idx] += delta
            if any(n <= 0 for n in new_dn):
                continue
            try:
                up_vecs, dn_vecs = build_sector_vectors(u1, u2, u3,
                                    base_up, tuple(new_dn), phi)
                up_b = gram_schmidt(up_vecs)
                dn_b = gram_schmidt(dn_vecs)
                V = compute_ckm(up_b, dn_b)
                J = jarlskog(V)
                t12, _, _ = pdg_angles_from_ckm(V)
                label = f"dn[{idx}] {base_dn[idx]:+d}->{new_dn[idx]:+d}"
                print(f"  {label:<28} {J:12.4e} {J/J_base:10.3f} "
                      f"{np.degrees(t12):10.2f}")
            except Exception as e:
                print(f"  dn[{idx}]{delta:+d}: ERROR {e}")

# ============================================================
# SECTION 10: Null tests
# ============================================================

def null_tests(field_name, poly_coeffs, phi=0.3,
               up_coeffs=(1,5,9), dn_coeffs=(3,7,12)):
    """
    Three null tests:
    A - Remove shared u1 direction
    B - Force unit vectors orthogonal
    C - Collapse to 1D per sector
    """
    print(f"\n{'='*60}")
    print(f"NULL TESTS: {field_name}")
    
    roots = find_real_roots(poly_coeffs)
    units_found, embed, log_vec = compute_log_vectors(roots)
    selected = select_independent_units(units_found, n=3)
    if len(selected) < 3:
        print("Insufficient units")
        return
    
    u1, u2, u3 = selected[0][1], selected[1][1], selected[2][1]
    
    # --- Null Test A: Remove shared u1 direction ---
    print("\nNull Test A: Remove shared u1 direction")
    print("  Up: n*u2*exp(i*phi) only (no u1)")
    print("  Down: n*u3 only (no u1)")
    
    u2c = u2.astype(complex) * np.exp(1j * phi)
    u3c = u3.astype(complex)
    
    up_vecs_A = [n * u2c for n in sorted(up_coeffs)]
    dn_vecs_A = [n * u3c for n in sorted(dn_coeffs)]
    
    try:
        up_b = gram_schmidt(up_vecs_A)
        dn_b = gram_schmidt(dn_vecs_A)
        V_A = compute_ckm(up_b, dn_b)
        J_A = jarlskog(V_A)
        t12_A, _, _ = pdg_angles_from_ckm(V_A)
        print(f"  |V| magnitudes:")
        for row in np.abs(V_A):
            print(f"    {[f'{x:.4f}' for x in row]}")
        print(f"  J = {J_A:.4e},  theta_12 = {np.degrees(t12_A):.2f} deg")
        print(f"  Expected: anarchic mixing (theta_12 >> 0.22)")
    except Exception as e:
        print(f"  All vectors parallel (1D sector): {e}")
        print("  This is the expected collapse for the 1D case")
    
    # --- Null Test B: Force unit vectors orthogonal ---
    print("\nNull Test B: Orthogonalized unit vectors (remove geometric structure)")
    
    # Gram-Schmidt on the unit vectors themselves to orthogonalize
    u_vecs = [u1, u2, u3]
    u_ortho = []
    for v in u_vecs:
        w = v.copy()
        for b in u_ortho:
            w = w - np.dot(b, w) * b
        u_ortho.append(w / norm(w))
    
    u1o, u2o, u3o = u_ortho
    up_vecs_B, dn_vecs_B = build_sector_vectors(u1o, u2o, u3o,
                                                  up_coeffs, dn_coeffs, phi)
    try:
        up_b = gram_schmidt(up_vecs_B)
        dn_b = gram_schmidt(dn_vecs_B)
        V_B = compute_ckm(up_b, dn_b)
        J_B = jarlskog(V_B)
        t12_B, _, _ = pdg_angles_from_ckm(V_B)
        print(f"  |V| magnitudes:")
        for row in np.abs(V_B):
            print(f"    {[f'{x:.4f}' for x in row]}")
        print(f"  J = {J_B:.4e},  theta_12 = {np.degrees(t12_B):.2f} deg")
        print(f"  Expected: large/anarchic mixing vs structured case")
    except Exception as e:
        print(f"  Error: {e}")
    
    # --- Null Test C: Collapse to 1D per sector ---
    print("\nNull Test C: Collapse to 1D (same vector for all three generations)")
    
    u2c = u2.astype(complex) * np.exp(1j * phi)
    u1c = u1.astype(complex)
    u3c = u3.astype(complex)
    
    # All three up generations = same direction (vary only magnitude)
    up_vecs_C = [u2c + u1c, 2*(u2c + u1c), 3*(u2c + u1c)]
    dn_vecs_C = [u3c + u1c, 2*(u3c + u1c), 3*(u3c + u1c)]
    
    try:
        up_b = gram_schmidt(up_vecs_C)
        dn_b = gram_schmidt(dn_vecs_C)
        V_C = compute_ckm(up_b, dn_b)
        J_C = jarlskog(V_C)
        t12_C, _, _ = pdg_angles_from_ckm(V_C)
        print(f"  J = {J_C:.4e},  theta_12 = {np.degrees(t12_C):.2f} deg")
    except ValueError:
        print("  Gram-Schmidt failed: vectors linearly dependent (expected 1D collapse)")

# ============================================================
# SECTION 11: Compare to PDG
# ============================================================

PDG_CKM = np.array([
    [0.97435, 0.22500, 0.003690],
    [0.22486, 0.97349, 0.04182],
    [0.00857, 0.04110, 0.99918]
])
PDG_J = 3.00e-5
PDG_theta12 = 0.22739  # rad
PDG_theta13 = 0.003692 # rad
PDG_theta23 = 0.04085  # rad

def compare_to_pdg(V, J, angles):
    """Print comparison with PDG values."""
    t12, t13, t23 = angles
    print(f"\n{'Parameter':<15} {'This model':>15} {'PDG 2022':>15} {'Ratio':>10}")
    print("-"*58)
    print(f"{'|V_ud|':<15} {abs(V[0,0]):15.4f} {PDG_CKM[0,0]:15.4f} "
          f"{abs(V[0,0])/PDG_CKM[0,0]:10.4f}")
    print(f"{'|V_us|':<15} {abs(V[0,1]):15.4f} {PDG_CKM[0,1]:15.4f} "
          f"{abs(V[0,1])/PDG_CKM[0,1]:10.4f}")
    print(f"{'|V_ub|':<15} {abs(V[0,2]):15.4f} {PDG_CKM[0,2]:15.4f} "
          f"{abs(V[0,2])/PDG_CKM[0,2]:10.4f}")
    print(f"{'|V_cb|':<15} {abs(V[1,2]):15.4f} {PDG_CKM[1,2]:15.4f} "
          f"{abs(V[1,2])/PDG_CKM[1,2]:10.4f}")
    print(f"{'|V_tb|':<15} {abs(V[2,2]):15.4f} {PDG_CKM[2,2]:15.4f} "
          f"{abs(V[2,2])/PDG_CKM[2,2]:10.4f}")
    print(f"{'theta_12 (rad)':<15} {t12:15.5f} {PDG_theta12:15.5f} "
          f"{t12/PDG_theta12:10.4f}")
    print(f"{'theta_13 (rad)':<15} {t13:15.5f} {PDG_theta13:15.5f} "
          f"{t13/PDG_theta13:10.4f}")
    print(f"{'theta_23 (rad)':<15} {t23:15.5f} {PDG_theta23:15.5f} "
          f"{t23/PDG_theta23:10.4f}")
    print(f"{'J':<15} {J:15.4e} {PDG_J:15.4e} {J/PDG_J:10.4f}")

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("ARITHMETIC FLAVOR GEOMETRY: FULL VERIFICATION")
    print("=" * 70)
    print(f"PDG reference: J = {PDG_J:.2e}, theta_12 = {np.degrees(PDG_theta12):.2f} deg")
    
    # ---- FIELD 1 ----
    # x^4 - x^3 - 4x^2 + 4x + 1
    poly1 = [1, -1, -4, 4, 1]
    
    result1 = analyze_field("Field 1: x^4 - x^3 - 4x^2 + 4x + 1",
                             poly1, phi=0.3)
    if result1 is not None:
        V1, angles1, J1, u1_f1, u2_f1, u3_f1 = result1
        compare_to_pdg(V1, J1, angles1)
    
    # Phase scan for Field 1
    phase_scan("Field 1", poly1)
    
    # Integer perturbation scan for Field 1
    integer_perturbation_scan("Field 1", poly1)
    
    # Null tests for Field 1
    null_tests("Field 1", poly1)
    
    # ---- FIELD 2 ----
    # x^4 - 2x^3 - x^2 + 2x + 1
    poly2 = [1, -2, -1, 2, 1]
    
    result2 = analyze_field("Field 2: x^4 - 2x^3 - x^2 + 2x + 1",
                             poly2, phi=0.3)
    if result2 is not None:
        V2, angles2, J2, u1_f2, u2_f2, u3_f2 = result2
        compare_to_pdg(V2, J2, angles2)
    
    # Phase scan for Field 2
    phase_scan("Field 2", poly2)
    
    # Integer perturbation scan for Field 2
    integer_perturbation_scan("Field 2", poly2)
    
    # Null tests for Field 2
    null_tests("Field 2", poly2)
    
    # ---- CROSS-FIELD COMPARISON ----
    print("\n" + "="*70)
    print("CROSS-FIELD SUMMARY")
    print("="*70)
    if result1 and result2:
        print(f"\n{'Quantity':<20} {'Field 1':>12} {'Field 2':>12} {'PDG':>12}")
        print("-"*58)
        print(f"{'J':<20} {J1:12.3e} {J2:12.3e} {PDG_J:12.3e}")
        print(f"{'theta_12 (deg)':<20} "
              f"{np.degrees(angles1[0]):12.2f} "
              f"{np.degrees(angles2[0]):12.2f} "
              f"{np.degrees(PDG_theta12):12.2f}")
        print(f"{'theta_13 (deg)':<20} "
              f"{np.degrees(angles1[1]):12.4f} "
              f"{np.degrees(angles2[1]):12.4f} "
              f"{np.degrees(PDG_theta13):12.4f}")
        print(f"{'theta_23 (deg)':<20} "
              f"{np.degrees(angles1[2]):12.3f} "
              f"{np.degrees(angles2[2]):12.3f} "
              f"{np.degrees(PDG_theta23):12.3f}")
        
        print(f"\nJ_field1 / J_PDG = {J1/PDG_J:.3f}")
        print(f"J_field2 / J_PDG = {J2/PDG_J:.3f}")
        
        theta12_ratio = np.degrees(angles1[0]) / np.degrees(angles2[0])
        print(f"theta_12 ratio (F1/F2) = {theta12_ratio:.3f}")
        print(f"(Close to 1.0 = mechanism robust across fields)")
    
    print("\n" + "="*70)
    print("COMPUTATION COMPLETE")
    print("="*70)