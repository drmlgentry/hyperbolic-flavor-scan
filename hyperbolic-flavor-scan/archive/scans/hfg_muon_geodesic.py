"""
Investigate the muon geodesic: ell ~ 5.293 in M_PMNS
q=44, ell_target = 44/4 * log(phi) = 11 * log(phi) = 5.29333

This geodesic matches m_mu/m_e = phi^11 to 0.003% precision.
What word in pi_1 realizes this geodesic?
What are its arithmetic properties?
"""
import snappy, numpy as np
from scipy.linalg import logm
import warnings; warnings.filterwarnings('ignore')

PHI = (1+5**0.5)/2
LOG_PHI = np.log(PHI)

M = snappy.OrientableClosedCensus[1]  # PMNS
rho = M.polished_holonomy()

print("="*65)
print("MUON GEODESIC INVESTIGATION")
print(f"Target: ell = 11*log(phi) = {11*LOG_PHI:.6f}")
print(f"Corresponds to: m_mu/m_e = phi^11 = {PHI**11:.4f} (PDG: 206.77)")
print("="*65)
print()

TARGET = 11 * LOG_PHI  # 5.29333
TOL = 0.01

def word_props(word):
    try:
        mat = np.array(rho(word), dtype=complex)
        mat = mat/np.sqrt(np.linalg.det(mat))
        ev = np.linalg.eigvals(mat)
        lam = ev[np.argmax(np.abs(ev))]
        ell = 2*abs(np.log(lam).real)
        twist = np.log(lam).imag
        h1 = (sum(1 for c in word if c=='a')-sum(1 for c in word if c=='A')+
              sum(1 for c in word if c=='b')-sum(1 for c in word if c=='B'))%5
        tr = float(np.real(np.trace(mat)))
        return ell, twist, h1, tr
    except: return None,None,None,None

def get_axis(word):
    try:
        mat = np.array(rho(word), dtype=complex)
        mat = mat/np.sqrt(np.linalg.det(mat))
        L = logm(mat)
        x = float(np.real(L[0,1]+L[1,0]))/2
        y = float(np.imag(L[1,0]-L[0,1]))/2
        z = float(np.real(L[0,0]-L[1,1]))/2
        v = np.array([x,y,z]); n = np.linalg.norm(v)
        return v/n if n>1e-10 else v
    except: return None

print("Searching words up to length 12 for geodesic near ell=5.293...")
print()

from itertools import product as iproduct
cancel = {'a':'A','A':'a','b':'B','B':'b'}

def valid(w):
    for i in range(len(w)-1):
        if cancel.get(w[i])==w[i+1]: return False
    return True

found = []
for length in range(1, 13):
    for combo in iproduct('abAB', repeat=length):
        w = ''.join(combo)
        if not valid(w): continue
        ell, twist, h1, tr = word_props(w)
        if ell is None: continue
        if abs(ell - TARGET) < TOL:
            found.append((abs(ell-TARGET), w, ell, twist, h1, tr))
    if length <= 8 or found:
        print(f"  Length {length}: {'found '+str(len(found))+' matches' if found else 'searching...'}")
    if found and length >= 8:
        break

print()
if found:
    found.sort()
    print(f"MATCHES FOUND (tol={TOL}):")
    print(f"  {'word':>15}  {'ell':>10}  {'dist':>8}  "
          f"{'twist_deg':>10}  {'H1':>4}  {'tr':>8}")
    print("  "+"-"*65)
    for dist, w, ell, twist, h1, tr in found[:10]:
        print(f"  {w:>15}  {ell:>10.6f}  {dist:>8.6f}  "
              f"{np.degrees(twist):>10.4f}  {h1:>4}  {tr:>8.4f}")
    
    # For the best match, compute detailed properties
    best = found[0]
    dist, w, ell, twist, h1, tr = best
    ax = get_axis(w)
    
    print()
    print(f"BEST MATCH: word='{w}'")
    print(f"  Geodesic length: {ell:.8f}")
    print(f"  Target (11*logphi): {TARGET:.8f}")
    print(f"  Distance: {dist:.2e}")
    print(f"  Relative error: {dist/TARGET*100:.6f}%")
    print(f"  Twist angle: {np.degrees(twist):.4f} deg")
    print(f"  H1 class: {h1} (mod 5)")
    print(f"  Trace: {tr:.6f}")
    print(f"  2*cosh(ell/2): {2*np.cosh(ell/2):.6f}")
    print(f"  phi^11+phi^-11 (L_11): {PHI**11+PHI**(-11):.6f}")
    print()
    print("  Lucas check: |tr| = L_11?")
    print(f"    |tr| = {abs(tr):.6f}")
    print(f"    L_11 = phi^11+phi^-11 = {PHI**11+PHI**(-11):.6f}")
    print(f"    2*cosh(ell/2) should = L_11 if ell=11*logphi exactly")
    print()
    
    # Check if ell = 11*log(phi) exactly via the bridge theorem
    print("  Bridge theorem check: ell = 2k*logphi <=> |tr| = L_k")
    print(f"    ell/(2*logphi) = {ell/(2*LOG_PHI):.6f}  (should be 11/2 = 5.5 if k=11/2)")
    print(f"    ell/logphi     = {ell/LOG_PHI:.6f}  (should be 11 if k=11)")
    print()
    
    # Axis alignment with PMNS mixing words
    print("  Axis alignment with PMNS mixing words:")
    mixing_words = ['aa', 'aaB', 'baa']
    ax_muon = get_axis(w)
    if ax_muon is not None:
        for mw in mixing_words:
            ax_mw = get_axis(mw)
            if ax_mw is not None:
                dot = float(np.dot(ax_muon, ax_mw))
                angle = np.degrees(np.arccos(np.clip(abs(dot), -1, 1)))
                print(f"    angle(muon_word, {mw}) = {angle:.4f} deg")

else:
    print("No matches found within tolerance.")
    print("Checking length spectrum directly...")
    ls = M.length_spectrum(cutoff=5.5)
    near = []
    for g in ls:
        try:
            v = g.length
            if callable(v): v=v()
            ell = float(complex(v).real)
            if abs(ell-TARGET) < 0.05:
                near.append((abs(ell-TARGET), ell))
        except: pass
    near.sort()
    print(f"Nearest geodesics to target {TARGET:.5f}:")
    for dist, ell in near[:5]:
        print(f"  ell={ell:.6f}  dist={dist:.6f}  "
              f"ell/logphi={ell/LOG_PHI:.4f}")

print()
print("="*65)
print("STRANGE QUARK GEODESIC CHECK")
print(f"Target: ell = 43/4 * log(phi) = {43/4*LOG_PHI:.6f}")
print("="*65)
print()

TARGET_S = 43/4 * LOG_PHI
ls = M.length_spectrum(cutoff=5.5)
near_s = []
for g in ls:
    try:
        v = g.length; 
        if callable(v): v=v()
        ell = float(complex(v).real)
        if abs(ell-TARGET_S) < 0.01:
            near_s.append((abs(ell-TARGET_S), ell))
    except: pass
near_s.sort()
if near_s:
    print(f"Matches for s-quark geodesic (q=43):")
    for dist, ell in near_s:
        print(f"  ell={ell:.6f}  dist={dist:.6f}  "
              f"ell/logphi={ell/LOG_PHI:.4f}  "
              f"relative error={dist/TARGET_S*100:.4f}%")
else:
    print("No match for s-quark geodesic within tol=0.01")
