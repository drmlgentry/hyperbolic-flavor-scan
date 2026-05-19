"""
slope_survey.py
===============
Map exp(2*ell_b - ell_a) across all Dehn fillings m003(p,q)
with |p|,|q| <= 8, and compare to known mass ratios.

Also check: does any linear combination c1*ell_a + c2*ell_b
with small integer coefficients give ln(known mass ratio)?

Run: conda run -n sage python slope_survey.py
"""
import snappy
import numpy as np

def get_ell(rho, word):
    """Real translation length from polished holonomy."""
    try:
        mat = np.array(rho(word), dtype=complex)
        mat /= np.sqrt(np.linalg.det(mat))
        eigs = np.linalg.eigvals(mat)
        lam = eigs[np.argmax(np.abs(eigs))]
        return 2*abs(float(np.real(np.log(lam))))
    except:
        return None

# Known mass ratios (PDG 2024)
mass_ratios = {
    'mb/mc':   4.18/1.27,
    'mc/ms':   1.27/0.095,
    'ms/md':   0.095/0.0047,
    'mt/mb':   172.69/4.18,
    'mb/ms':   4.18/0.095,
    'mtau/mmu': 1.7769/0.10566,
    'mmu/me':  0.10566/0.000511,
}

print("="*65)
print("SLOPE SURVEY: exp(c1*ell_a + c2*ell_b) on m003(p,q)")
print("="*65)
print()
print("Checking combinations (c1,c2) with |c1|,|c2| <= 3")
print()

m003_cusped = snappy.Manifold("m003")

# First establish the target slope (-2,3) baseline
print("TARGET SLOPE (-2,3) BASELINE:")
print("-"*50)
M_target = snappy.OrientableClosedCensus[1]
rho_target = M_target.polished_holonomy()
ell_a_t = get_ell(rho_target, 'a')
ell_b_t = get_ell(rho_target, 'b')
print(f"  ell_a = {ell_a_t:.8f}")
print(f"  ell_b = {ell_b_t:.8f}")
print()

# Check all linear combinations at target slope
print("Linear combinations at slope (-2,3):")
for c1 in range(-3, 4):
    for c2 in range(-3, 4):
        if c1 == 0 and c2 == 0: continue
        val = c1*ell_a_t + c2*ell_b_t
        if val <= 0: continue
        r = np.exp(val)
        for name, target in mass_ratios.items():
            err = abs(r - target)/target*100
            if err < 2.0:
                print(f"  ({c1:>2},{c2:>2}): exp({c1}*ell_a + {c2}*ell_b)"
                      f" = exp({val:.6f}) = {r:.6f}"
                      f" ~ {name}={target:.6f}  err={err:.4f}%")

print()
print("="*65)
print("SLOPE SURVEY: mapping exp(2*ell_b - ell_a) across slopes")
print("="*65)
print()
print(f"{'slope':>10}  {'ell_a':>10}  {'ell_b':>10}  "
      f"{'2b-a':>10}  {'exp(2b-a)':>10}  {'err%':>8}  note")
print("-"*80)

MB_MC = 4.18/1.27
hits = []

for p in range(-8, 9):
    for q in range(1, 9):  # positive q by convention
        if p == 0: continue
        try:
            M = m003_cusped.copy()
            M.dehn_fill((p, q))
            if not 'positively' in str(M.solution_type()):
                continue
            rho = M.polished_holonomy()
            ell_a = get_ell(rho, 'a')
            ell_b = get_ell(rho, 'b')
            if ell_a is None or ell_b is None: continue
            if ell_a < 0.01 or ell_b < 0.01: continue
            combo = 2*ell_b - ell_a
            if combo <= 0: continue
            r = np.exp(combo)
            err = abs(r - MB_MC)/MB_MC*100
            flag = ""
            if p == -2 and q == 3: flag = " <-- TARGET"
            if err < 5.0:
                hits.append((err, p, q, ell_a, ell_b, combo, r))
            print(f"  ({p:>3},{q:>2}):  {ell_a:>10.6f}  {ell_b:>10.6f}  "
                  f"{combo:>10.6f}  {r:>10.6f}  {err:>8.3f}%{flag}")
        except:
            pass

print()
print("="*65)
print("SLOPES WHERE exp(2*ell_b-ell_a) IS WITHIN 5% OF mb/mc")
print("="*65)
hits.sort()
for err, p, q, ea, eb, combo, r in hits:
    flag = " <-- TARGET" if p==-2 and q==3 else ""
    print(f"  ({p:>3},{q:>2}): err={err:.4f}%  "
          f"ell_a={ea:.5f}  ell_b={eb:.5f}  "
          f"exp={r:.6f}{flag}")

print()
print("="*65)
print("m006 SLOPE SURVEY: same analysis")
print("="*65)
print()

m006_cusped = snappy.Manifold("m006")
M_ckm = snappy.OrientableClosedCensus[43]
rho_ckm = M_ckm.polished_holonomy()
ell_a_c = get_ell(rho_ckm, 'a')
ell_b_c = get_ell(rho_ckm, 'b')
print(f"m006(-5,2) baseline: ell_a={ell_a_c:.8f}, ell_b={ell_b_c:.8f}")
print()
print("Linear combinations at slope (-5,2):")
for c1 in range(-3, 4):
    for c2 in range(-3, 4):
        if c1 == 0 and c2 == 0: continue
        val = c1*ell_a_c + c2*ell_b_c
        if val <= 0: continue
        r = np.exp(val)
        for name, target in mass_ratios.items():
            err = abs(r - target)/target*100
            if err < 2.0:
                print(f"  ({c1:>2},{c2:>2}): exp({c1}*ell_a + {c2}*ell_b)"
                      f" = {r:.6f} ~ {name}={target:.6f}  err={err:.4f}%")
