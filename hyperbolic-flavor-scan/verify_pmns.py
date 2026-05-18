"""
verify_pmns.py
Verify all canonical PMNS values needed for the paper rebuild.
Run in WSL: conda run -n sage python verify_pmns.py
"""
import snappy
import numpy as np
from scipy.linalg import logm
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("PMNS CANONICAL VALUE VERIFICATION")
print("="*60)

# M_PMNS must be loaded from WSL sage environment
M = snappy.OrientableClosedCensus[1]
print(f"\nManifold: {M.name()}")
print(f"Volume:   {M.volume():.5f}")
print(f"H1:       {M.homology()}")
try:
    cs = M.chern_simons()
    print(f"CS:       {cs:.5f}")
except:
    print(f"CS:       (use dehn_surgery_equations or literature: 1/4)")

rho = M.polished_holonomy()

# ── H1 classes ─────────────────────────────────────────────────────
print("\nH1 classes ([a]=?, [b]=? in H1=Z/5):")
print("  (Reading from abelianization map)")
# The homology map: compute [word] mod 5
# Try to determine from the manifold's peripheral structure
try:
    ab = M.abelian_cover(2)
    print(f"  degree-2 abelian cover: {ab}")
except:
    pass

# Compute from word traces (empirical check)
def word_trace(word):
    mat = np.array(rho(word), dtype=complex)
    return np.trace(mat)

print()
print("Holonomy traces of key words:")
for word in ['a','b','aa','ab','ba','aaB','baa','aab','bba']:
    try:
        tr = word_trace(word)
        print(f"  tr(rho({word})) = {tr.real:.6f} + {tr.imag:.6f}i")
    except:
        print(f"  {word}: failed")

# ── Axis extraction ─────────────────────────────────────────────────
print()
print("Axis vectors for canonical triple {aa, aaB, baa}:")

def get_axis_and_twist(word):
    mat = np.array(rho(word), dtype=complex)
    mat /= np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    nx = float(np.real(L[0,1]+L[1,0]))/2
    ny = float(np.imag(L[1,0]-L[0,1]))/2
    nz = float(np.real(L[0,0]-L[1,1]))/2
    v = np.array([nx,ny,nz])
    norm = np.linalg.norm(v)
    axis = v/norm if norm > 1e-10 else v
    # Twist angle: Im(log(lambda)) where lambda = dominant eigenvalue
    evals = np.linalg.eigvals(mat)
    lam = evals[0] if abs(evals[0]) >= abs(evals[1]) else evals[1]
    twist = float(np.imag(np.log(lam)))
    return axis, twist

axes = {}
twists = {}
for word in ['aa','aaB','baa']:
    try:
        ax, tw = get_axis_and_twist(word)
        axes[word] = ax
        twists[word] = tw
        print(f"  {word}: n={np.round(ax,5)}, phi={np.degrees(tw):.4f} deg")
    except Exception as e:
        print(f"  {word}: FAILED {e}")

# Pairwise angles
print()
print("Pairwise axis angles:")
words = ['aa','aaB','baa']
for i in range(3):
    for j in range(i+1,3):
        w1,w2 = words[i],words[j]
        if w1 in axes and w2 in axes:
            ang = np.degrees(np.arccos(
                np.clip(abs(np.dot(axes[w1],axes[w2])),-1,1)))
            print(f"  {w1}-{w2}: {ang:.4f} deg")

# ── CP phase ────────────────────────────────────────────────────────
print()
print("CP phase formula: delta = pi + phi(aaB) + phi(baa)")
if 'aaB' in twists and 'baa' in twists:
    phi_aaB = np.degrees(twists['aaB'])
    phi_baa = np.degrees(twists['baa'])
    delta = 180 + phi_aaB + phi_baa
    print(f"  phi(aaB) = {phi_aaB:.4f} deg")
    print(f"  phi(baa) = {phi_baa:.4f} deg")
    print(f"  delta    = {delta:.4f} deg")
    print(f"  PDG:       197.0 deg")
    print(f"  Error:     {abs(delta-197.0)/197.0*100:.4f}%")

# ── Sigma opt geodesic ──────────────────────────────────────────────
print()
phi_gold = (1+5**0.5)/2
sigma_geod = 1.5*np.log(phi_gold)
print(f"Geodesic claim: sigma_opt = 3/2*log(phi) = {sigma_geod:.5f}")
print("Checking geodesic length spectrum:")
try:
    spec = M.length_spectrum(1.5)
    for g in spec[:8]:
        print(f"  length={g.length.real:.5f}, "
              f"ratio to log(phi)={g.length.real/np.log(phi_gold):.4f}")
except Exception as e:
    print(f"  length_spectrum failed: {e}")

# ── Lepton norm witnesses ───────────────────────────────────────────
print()
print("Lepton Z[omega] norm witnesses:")
def eisenstein_norm(a,b):
    return a*a - a*b + b*b

for name,(ratio,a,b,claimed) in [
    ('mu',  (206.77,  -16, -12, 208)),
    ('tau', (3477.22, -68, -37, 3477)),
]:
    N = eisenstein_norm(a,b)
    err = abs(ratio-N)/ratio*100
    ok = N==claimed
    print(f"  {name}: N({a},{b})=a^2-ab+b^2={N} [claim={claimed}] "
          f"{'OK' if ok else 'WRONG'}, err={err:.4f}%")

print()
print("="*60)
print("Run complete.")
