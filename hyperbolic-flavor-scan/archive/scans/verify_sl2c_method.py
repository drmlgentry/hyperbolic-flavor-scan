"""
verify_sl2c_method.py
Re-verify twist paper claims using the ORIGINAL method from twist_census.py:
G.SL2C(word) instead of polished_holonomy()
Also verify the phi_fold convention and MZ/MW.
Run: conda run -n sage python verify_sl2c_method.py
"""
import snappy, numpy as np

M_P = snappy.OrientableClosedCensus[1]
M_C = snappy.OrientableClosedCensus[43]

def get_eigendata_sl2c(G, word):
    """Original method from twist_census.py using G.SL2C()."""
    try:
        mat = np.array([[complex(G.SL2C(word)[i,j]) for j in range(2)]
                        for i in range(2)])
        eigs = np.linalg.eigvals(mat)
        # Original convention: take eigenvalue with Im(log) >= 0
        lam = eigs[0] if np.imag(eigs[0]) >= 0 else eigs[1]
        log_lam = np.log(lam)
        t = np.real(log_lam)
        phi = np.imag(log_lam)  # radians
        return t, phi, np.exp(abs(t)), np.degrees(phi)
    except:
        return None, None, None, None

def phi_fold(phi_deg):
    """Various folding conventions used in original script."""
    return {
        'direct': phi_deg,
        'abs': abs(phi_deg),
        '180-abs': 180 - abs(phi_deg),
        'mod180': phi_deg % 180,
        'abs_mod90': abs(phi_deg) % 90,
        '360-abs': 360 - abs(phi_deg),
    }

from itertools import product as iprod

print("="*60)
print("VERIFICATION USING ORIGINAL SL2C METHOD")
print("="*60)
print()

# SM targets
sm_targets = {
    'theta13_NU': 8.54,
    'theta13_CKM': 0.201,
    'MZ_MW': 1.13451 * 90,  # as an angle? No -- check mass ratio
    'delta_CKM': 68.0,
    'theta12_NU': 33.41,
}

for manifold, M, name in [
    ('m003 (PMNS)', M_P, 'm003'),
    ('m006 (CKM)',  M_C, 'm006'),
]:
    print(f"\n{'='*50}")
    print(f"{manifold}")
    print(f"{'='*50}")
    G = M.fundamental_group()
    ngen = len(list(G.generators()))
    letters = [c for c in 'abcd'[:ngen]] + [c.upper() for c in 'abcd'[:ngen]]

    print(f"Generators: {list(G.generators())}")
    print()

    # Scan words up to length 4
    print("Searching for theta13_NU=8.54 and theta13_CKM=0.201 (within 0.5 deg):")
    hits_t13 = []
    for length in range(1, 5):
        for combo in iprod(letters, repeat=length):
            w = ''.join(combo)
            t, phi_rad, mod_lam, phi_deg = get_eigendata_sl2c(G, w)
            if t is None: continue
            # Try all folding conventions
            for fold_name, val in phi_fold(phi_deg).items():
                for target_name, target in [('theta13_NU', 8.54),
                                             ('theta13_CKM', 0.201)]:
                    if abs(val - target) < 0.5:
                        hits_t13.append((abs(val-target), w, fold_name,
                                         val, phi_deg, target_name))

    hits_t13.sort()
    for err, w, fold, val, phi_raw, tname in hits_t13[:10]:
        print(f"  {tname}: {w} [{fold}({phi_raw:.3f}°)] = {val:.4f}° err={err:.4f}°")

    # Also look for MZ/MW as ratio of mod_lambda
    print()
    print("Searching for MZ/MW=1.1345 as |lambda| ratio (within 1%):")
    lams = {}
    for length in range(1, 5):
        for combo in iprod(letters, repeat=length):
            w = ''.join(combo)
            t, phi_rad, mod_lam, phi_deg = get_eigendata_sl2c(G, w)
            if t is not None and mod_lam > 1.0:
                lams[w] = (mod_lam, abs(t)*2, phi_deg)

    MZ_MW = 1.13451
    hits_mzmw = []
    wlist = list(lams.items())
    for i, (w1, (l1, e1, p1)) in enumerate(wlist):
        for w2, (l2, e2, p2) in wlist[i+1:]:
            for r, meth in [(l1/l2, 'lam'), (e1/e2, 'ell')]:
                if abs(r - MZ_MW) < 0.02:
                    hits_mzmw.append((abs(r-MZ_MW), w1, w2, r, meth))
    hits_mzmw.sort()
    for err, w1, w2, r, meth in hits_mzmw[:5]:
        print(f"  {meth}: {w1}/{w2} = {r:.6f} err={err:.6f}")

    # Verify known claims
    print()
    print("Verifying known claims:")
    known = {
        'm003': [
            ('aa',  'CP phi(aa)',          168.5, 0.5),
            ('aaB', 'CP phi(aaB)',         -176.7, 0.5),
            ('baa', 'CP phi(baa)',         -167.4, 0.5),
            ('bbbb', 'mb/mc numerator',    None, None),
            ('bAbA', 'mb/mc denominator',  None, None),
            ('AAB', 'theta12_CKM fold',    12.6, 0.5),
        ],
        'm006': [
            ('aa',   'delta_CKM fold',     67.6, 0.5),
            ('aaB',  'isospectral1',       -87.5, 0.5),
            ('AbA',  'isospectral2',       -87.5, 0.5),
            ('AAb',  'isospectral3',       -87.5, 0.5),
            ('aaabb','theta23_CKM',         2.13, 0.5),
            ('abbAB','solar_angle',        33.6, 0.5),
        ],
    }
    for w, desc, expected, tol in known.get(name, []):
        t, phi_rad, mod_lam, phi_deg = get_eigendata_sl2c(G, w)
        if t is None:
            print(f"  {w} ({desc}): FAILED to compute")
            continue
        if expected is None:
            print(f"  {w} ({desc}): phi={phi_deg:.4f}°, |lam|={mod_lam:.6f}")
        else:
            match = abs(phi_deg - expected) < tol
            status = "OK" if match else f"OFF by {abs(phi_deg-expected):.3f}"
            print(f"  {w} ({desc}): phi={phi_deg:.4f}° [{status}]")

print()
print("="*60)
print("mb/mc ratio check (SL2C method):")
print("="*60)
G_P = M_P.fundamental_group()
results = {}
for w in ['bbbb', 'bAbA']:
    t, phi_rad, mod_lam, phi_deg = get_eigendata_sl2c(G_P, w)
    results[w] = mod_lam
    print(f"  |lambda({w})| = {mod_lam:.6f}")
if results.get('bbbb') and results.get('bAbA'):
    ratio = results['bbbb']/results['bAbA']
    mb_mc = 4.18/1.27
    print(f"  Ratio = {ratio:.6f}")
    print(f"  PDG mb/mc = {mb_mc:.6f}")
    print(f"  Error = {abs(ratio-mb_mc)/mb_mc*100:.4f}%")
