import numpy as np
from itertools import product

def find_totally_real_cubic_roots(coeffs):
    """Find roots of x^3 + ax^2 + bx + c, check totally real."""
    roots = np.roots([1] + coeffs)
    if np.all(np.abs(roots.imag) < 1e-10):
        return sorted(roots.real)
    return None

def log_embedding(unit_val, embeddings):
    """Compute log embedding vector for a unit given its values at embeddings."""
    return np.array([np.log(abs(e)) for e in embeddings])

def find_units(roots, max_coeff=8):
    """
    Search for units in Z[alpha] by scanning a+b*alpha+c*alpha^2.
    A unit has norm +-1 (product of all conjugate values = +-1).
    """
    r1, r2, r3 = roots
    units = []
    seen = set()
    for a, b, c in product(range(-max_coeff, max_coeff+1), repeat=3):
        if a == 0 and b == 0 and c == 0:
            continue
        # Value at each embedding
        v1 = a + b*r1 + c*r1**2
        v2 = a + b*r2 + c*r2**2
        v3 = a + b*r3 + c*r3**2
        norm = v1 * v2 * v3
        if abs(abs(norm) - 1.0) < 1e-6:
            if abs(v1) > 1e-9 and abs(v2) > 1e-9 and abs(v3) > 1e-9:
                # Canonical form: use log embedding as key
                lv = tuple(round(x, 6) for x in sorted([
                    np.log(abs(v1)), np.log(abs(v2)), np.log(abs(v3))
                ]))
                if lv not in seen and not all(abs(x) < 1e-8 for x in lv):
                    seen.add(lv)
                    units.append((a, b, c, v1, v2, v3))
    return units

def gram_matrix_and_kappa(units, roots):
    """
    Pick two independent unit log vectors, compute Gram matrix,
    return eigenvalues and anisotropy ratio kappa.
    """
    vecs = []
    for (a, b, c, v1, v2, v3) in units:
        lv = np.array([np.log(abs(v1)), np.log(abs(v2)), np.log(abs(v3))])
        vecs.append(lv)

    # Find two linearly independent vectors
    basis = []
    for v in vecs:
        if len(basis) == 0:
            if np.linalg.norm(v) > 1e-8:
                basis.append(v)
        elif len(basis) == 1:
            # Check independence
            cross = np.cross(basis[0][:2] if len(basis[0]) > 2 else basis[0],
                           v[:2] if len(v) > 2 else v)
            mat = np.array([basis[0], v])
            if np.linalg.matrix_rank(mat, tol=1e-6) == 2:
                basis.append(v)
                break

    if len(basis) < 2:
        return None, None, None

    # Gram matrix (project to sum=0 hyperplane implicitly via dot products)
    G = np.array([[np.dot(basis[i], basis[j]) for j in range(2)]
                  for i in range(2)])
    eigvals = np.sort(np.linalg.eigvalsh(G))
    eigvals = eigvals[eigvals > 1e-10]
    if len(eigvals) < 2:
        return G, eigvals, None
    kappa = eigvals[-1] / eigvals[0]
    return G, eigvals, kappa

# ============================================================
# Test fields: totally real cubics
# Polynomial x^3 + ax^2 + bx + c given as [a, b, c]
# ============================================================
fields = [
    ([0, -3, 1],   "x^3 - 3x + 1          (disc=81)"),
    ([-1, -2, 1],  "x^3 - x^2 - 2x + 1    (disc=49)"),
    ([0, -4, 1],   "x^3 - 4x + 1          (disc=229)"),
    ([-1, -3, 1],  "x^3 - x^2 - 3x + 1    (disc=148)"),
    ([0, -5, 1],   "x^3 - 5x + 1          (disc=473)"),
    ([-2, -1, 1],  "x^3 - 2x^2 - x + 1    (disc=49)"),
    ([0, -6, 2],   "x^3 - 6x + 2          (disc=*)"),
    ([-1, -4, 3],  "x^3 - x^2 - 4x + 3    (disc=*)"),
    ([-3, -1, 1],  "x^3 - 3x^2 - x + 1    (disc=*)"),
    ([-1, -5, 2],  "x^3 - x^2 - 5x + 2    (disc=*)"),
]

print("=" * 70)
print("REGULATOR LATTICE ANISOTROPY: TOTALLY REAL CUBIC FIELDS")
print("=" * 70)
print(f"{'Field':<40} {'Units':>6} {'lam1':>8} {'lam2':>8} {'kappa':>8}")
print("-" * 70)

kappas = []
for coeffs, label in fields:
    roots = find_totally_real_cubic_roots(coeffs)
    if roots is None:
        print(f"{label:<40} {'NOT TOTALLY REAL':>30}")
        continue
    units = find_units(roots, max_coeff=6)
    if len(units) < 2:
        print(f"{label:<40} {'<2 units found':>30}")
        continue
    G, eigvals, kappa = gram_matrix_and_kappa(units, roots)
    if kappa is None:
        print(f"{label:<40} {len(units):>6} {'indep fail':>18}")
        continue
    kappas.append(kappa)
    print(f"{label:<40} {len(units):>6} {eigvals[0]:>8.3f} {eigvals[-1]:>8.3f} {kappa:>8.3f}")

print("-" * 70)
if kappas:
    print(f"kappa range: {min(kappas):.3f} -- {max(kappas):.3f}")
    print(f"mean kappa:  {np.mean(kappas):.3f}")
    print(f"max kappa:   {max(kappas):.3f}")
    print()
    if max(kappas) < 20:
        print("RESULT: kappa bounded across all tested fields.")
        print("Supports Conjecture 1 (bounded anisotropy).")
    else:
        print("RESULT: large kappa found -- conjecture may need revision.")

print()
print("Spectrum obstruction reminder:")
print("  Need kappa >= 34.46/N for N-bounded coordinates.")
print("  For N=6: kappa >= 5.74 required to fit 9 fermion masses.")
print("  If max kappa << 5.74, cubic regulator lattice is ruled out.")
