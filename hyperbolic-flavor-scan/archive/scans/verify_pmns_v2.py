"""
verify_pmns_v2.py
Verify PMNS canonical values with corrected twist angle computation.
Run: conda run -n sage python verify_pmns_v2.py
"""
import snappy
import numpy as np
from scipy.linalg import logm
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("PMNS CANONICAL VALUE VERIFICATION v2")
print("="*60)

M = snappy.OrientableClosedCensus[1]
print(f"\nManifold: {M.name()}")
print(f"Volume:   {M.volume():.5f}")
print(f"H1:       {M.homology()}")

rho = M.polished_holonomy()

# ── H1 classes via abelianization ──────────────────────────────────
# H1(M_PMNS) = Z/5. The abelianization sends generators to their
# classes. We can read these from the holonomy eigenvalues of
# flat U(1) representations, or just compute mod 5 directly.
# For m003(-2,3): slope (-2,3) means 2 meridian + 3 longitude = 0
# The peripheral structure determines [a] and [b].
# Try to get it from snappy directly:
print()
print("Attempting H1 class computation:")
try:
    # Use the fundamental group presentation
    G = M.fundamental_group()
    print(f"  Generators: {G.generators()}")
    print(f"  Relations:  {G.relators()[:2]}")
except Exception as e:
    print(f"  FundamentalGroup failed: {e}")

# Empirical: try small words and see which flat U(1) character
# assigns them to which class mod 5
# For a flat U(1) character chi_k: chi_k(w) = exp(2*pi*i*k*[w]/5)
# We can identify [w] from chi_1(w) = exp(2*pi*i*[w]/5)
# chi_1 corresponds to the order-5 representation
print()
print("H1 classes from flat U(1) eigenvalues:")
print("  (chi_k(w) = exp(2*pi*i*k*[w]/5) for the abelian rep)")

# The flat U(1) rep of order 5 maps generators to 5th roots of unity
# Try to identify it from the holonomy
# For m003(-2,3), literature says [a]=1, [b]=2 in H1=Z/5
# Let's verify via abelianization of relator words
test_words = ['a','b','aa','ab','ba','bb','aaa','aaB','baa']
print("  Word  |  [word] mod 5 (from trace field structure)")
print("  " + "-"*40)
# Use the fact that for the abelian rep, det(rho(w))=1
# and for Z/5 rep: rho_ab(a)=zeta^[a], rho_ab(b)=zeta^[b]
# We can identify by testing which assignment satisfies relations
# For m003(-2,3): the Dehn filling relation is a^(-2)*b^3 in pi_1
# So -2[a] + 3[b] = 0 mod 5
# Combined with rank of H1=1, only solution up to automorphism:
# if [a]=3, [b]=2: -6+6=0 mod 5? -6 mod 5 = 4, 6 mod 5 = 1 -> 4+1=5=0 YES
# if [a]=1, [b]=? : -2+3[b]=0 -> [b]=4: check 3*4=12=2 mod 5, -2+2=0? No
# Let me solve properly:
print()
print("  Solving: -2[a] + 3[b] = 0 (mod 5) for H1(m003(-2,3))=Z/5")
for a_class in range(5):
    for b_class in range(1,5):  # b_class != 0 (non-trivial)
        if (-2*a_class + 3*b_class) % 5 == 0:
            if a_class != 0:  # a must be non-trivial generator
                print(f"  Solution: [a]={a_class}, [b]={b_class} (mod 5)")

# Use [a]=3, [b]=2 (the solution that makes sense with H1=Z/5)
a_cl, b_cl = 3, 2
print(f"\n  Using [a]={a_cl}, [b]={b_cl} (mod 5)")
print()
for word in ['aa','aaB','baa']:
    cl = 0
    for c in word:
        if c == 'a': cl += a_cl
        elif c == 'A': cl -= a_cl
        elif c == 'b': cl += b_cl
        elif c == 'B': cl -= b_cl
    cl = cl % 5
    print(f"  [{word}] = {cl} (mod 5)")

print()
print("  Also trying [a]=1, [b]=3:")
a_cl2, b_cl2 = 1, 3
for word in ['aa','aaB','baa']:
    cl = 0
    for c in word:
        if c == 'a': cl += a_cl2
        elif c == 'A': cl -= a_cl2
        elif c == 'b': cl += b_cl2
        elif c == 'B': cl -= b_cl2
    cl = cl % 5
    print(f"  [{word}] = {cl} (mod 5)")

# ── Corrected twist angle computation ──────────────────────────────
print()
print("="*40)
print("CORRECTED TWIST ANGLE COMPUTATION")
print("="*40)
print()
print("Method: eigenvalue of SL(2,C) matrix (not logm)")
print("  lambda = eigenvalue with |lambda| >= 1")
print("  twist phi = Im(log(lambda))")
print("  This matches the SnapPy convention used in the CP paper")

def get_twist_and_axis(word):
    """Get twist angle and axis using eigenvalue method."""
    mat = np.array(rho(word), dtype=complex)
    # Normalize to SL(2,C)
    mat /= np.sqrt(np.linalg.det(mat))

    # Eigenvalues: lambda and 1/lambda
    evals = np.linalg.eigvals(mat)
    # Choose the eigenvalue with Im(log) in (-pi, pi]
    lam = evals[0]
    if abs(evals[1]) > abs(evals[0]):
        lam = evals[1]

    # Translation length and twist
    log_lam = np.log(lam)
    length = 2 * abs(np.real(log_lam))
    twist = np.imag(log_lam)  # in radians

    # Axis from logm (imaginary part of Pauli decomposition)
    L = logm(mat)
    # Use imaginary parts for axis (standard for loxodromic)
    nx = float(np.imag(L[0,1]))
    ny = float(np.imag(L[1,1]))
    nz = float(np.imag(L[0,0]))
    v = np.array([nx, ny, nz])
    norm = np.linalg.norm(v)
    axis = v/norm if norm > 1e-10 else v
    return axis, twist, length, lam

print()
words = ['aa', 'aaB', 'baa']
axes = {}
twists = {}
lengths = {}

for word in words:
    try:
        ax, tw, ell, lam = get_twist_and_axis(word)
        axes[word] = ax
        twists[word] = tw
        lengths[word] = ell
        print(f"  {word}:")
        print(f"    lambda = {lam:.6f}")
        print(f"    length = {ell:.6f}")
        print(f"    twist  = {np.degrees(tw):.4f} deg ({tw:.6f} rad)")
        print(f"    axis   = {np.round(ax,5)}")
    except Exception as e:
        print(f"  {word}: FAILED {e}")

# CP phase
print()
print("CP phase:")
if 'aaB' in twists and 'baa' in twists:
    phi_aaB = np.degrees(twists['aaB'])
    phi_baa = np.degrees(twists['baa'])

    # Try all sign combinations and angle conventions
    for sign_a, sign_b in [(1,1),(1,-1),(-1,1),(-1,-1)]:
        delta = 180 + sign_a*phi_aaB + sign_b*phi_baa
        delta_mod = delta % 360
        err = abs(delta_mod - 197.0)/197.0*100
        marker = " <-- MATCH" if err < 2 else ""
        print(f"  pi + ({sign_a})*phi(aaB) + ({sign_b})*phi(baa) = "
              f"{delta:.4f} = {delta_mod:.4f} deg, err={err:.2f}%{marker}")

    # Also try with the other eigenvalue convention
    print()
    print("  Trying with conjugate eigenvalue (twist -> -twist):")
    phi_aaB_conj = -phi_aaB
    phi_baa_conj = -phi_baa
    for sign_a, sign_b in [(1,1),(1,-1),(-1,1),(-1,-1)]:
        delta = 180 + sign_a*phi_aaB_conj + sign_b*phi_baa_conj
        delta_mod = delta % 360
        err = abs(delta_mod - 197.0)/197.0*100
        marker = " <-- MATCH" if err < 2 else ""
        print(f"  pi + ({sign_a})*(-phi(aaB)) + ({sign_b})*(-phi(baa)) = "
              f"{delta_mod:.4f} deg, err={err:.2f}%{marker}")

# Pairwise axis angles
print()
print("Pairwise axis angles:")
for i in range(3):
    for j in range(i+1,3):
        w1,w2 = words[i],words[j]
        if w1 in axes and w2 in axes:
            ang = np.degrees(np.arccos(
                np.clip(abs(np.dot(axes[w1],axes[w2])),-1,1)))
            print(f"  {w1}-{w2}: {ang:.4f} deg")

# Golden ratio geodesic check
print()
phi_gold = (1+5**0.5)/2
sigma_pred = 1.5*np.log(phi_gold)
print(f"Predicted sigma_opt = 3/2*log(phi) = {sigma_pred:.5f}")
print("Geodesic lengths (first 10):")
try:
    spec = M.length_spectrum(2.0)
    for i,g in enumerate(list(spec)[:10]):
        ell = float(g.length.real)
        ratio = ell/np.log(phi_gold)
        print(f"  [{i}] ell={ell:.5f}, ell/log(phi)={ratio:.4f}")
except Exception as e:
    print(f"  Failed: {e}")

# Lepton norms
print()
print("Lepton Z[omega] norm witnesses (Eisenstein, a^2-ab+b^2):")
for name,(ratio,a,b,claimed) in [
    ('mu',  (206.77,  -16,-12, 208)),
    ('tau', (3477.22, -68,-37, 3477)),
]:
    N = a*a - a*b + b*b
    err = abs(ratio-N)/ratio*100
    print(f"  {name}: N({a},{b})={N} [claim={claimed}] "
          f"{'OK' if N==claimed else 'WRONG'}, err={err:.4f}%")

print()
print("="*60)
