"""
Fix: verify k=11 (L_11=199) is Eisenstein norm
Identify CKM trace field via the cusped ancestor m006
"""
import snappy, numpy as np
import warnings; warnings.filterwarnings('ignore')

PHI = (1+5**0.5)/2
LOG_PHI = np.log(PHI)
m_e = 0.000511

def is_prime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5)+1, 2):
        if n%i==0: return False
    return True

def find_all_eis_norms(N, limit=300):
    """Brute force: find all (a,b) with a^2-ab+b^2=N"""
    results = []
    for a in range(-limit, limit+1):
        for b in range(-limit, limit+1):
            if a*a - a*b + b*b == N:
                results.append((a,b))
    return results

def find_all_gauss_norms(N, limit=300):
    """Brute force: find all (a,b) with a^2+b^2=N"""
    results = []
    for a in range(0, int(N**0.5)+1):
        b_sq = N - a*a
        b = int(b_sq**0.5+0.5)
        if b*b == b_sq:
            results.append((a,b))
            if a != b and b != 0:
                results.append((b,a))
    return results

# ── Fix: brute-force check for k=11 ─────────────────────────────────────────
print("="*65)
print("FIX: BRUTE-FORCE EISENSTEIN NORM CHECK FOR k=11")
print("="*65)
print()

for k in [4, 11, 17, 19]:
    Lk = PHI**k + PHI**(-k)
    Lk_int = round(Lk)
    eis = find_all_eis_norms(Lk_int, limit=300)
    gauss = find_all_gauss_norms(Lk_int)
    mass = m_e * PHI**k

    print(f"k={k}: L_k={Lk:.4f}, L_k_int={Lk_int}")
    print(f"  prime: {is_prime(Lk_int)}")
    print(f"  mod 3: {Lk_int%3} ({'splits' if Lk_int%3==1 else 'inert/ram'})")
    print(f"  Eisenstein norms (a^2-ab+b^2={Lk_int}): "
          f"{len(eis)} solutions, e.g. {eis[:3]}")
    print(f"  Gaussian norms (a^2+b^2={Lk_int}): "
          f"{len(gauss)} solutions, e.g. {gauss[:3]}")
    print(f"  phi^k = {PHI**k:.6f}")
    print(f"  mass  = m_e * phi^k = {mass:.4f} GeV")

    # SM match
    sm = {'mu':0.10566,'tau':1.77686,'c':1.275,'b':4.18}
    for name, fmass in sm.items():
        err = abs(mass-fmass)/fmass*100
        if err < 20:
            print(f"  SM match: {name} = {fmass} GeV ({err:.2f}% error)")
    print()

# ── Full corrected table ─────────────────────────────────────────────────────
print("="*65)
print("CORRECTED LUCAS-EISENSTEIN TABLE (brute force)")
print("="*65)
print()
print(f"{'k':>4} {'L_k_int':>10} {'prime':>6} {'mod3':>5} "
      f"{'EisNorm':>8} {'GaussNorm':>10} {'mass_GeV':>10} {'SM':>15}")
print("-"*75)

all_results = []
for k in range(0, 40):
    Lk = PHI**k + PHI**(-k)
    Lk_int = round(Lk)
    if abs(Lk - Lk_int) > 0.015:
        continue
    if Lk_int < 2:
        continue
    prime = is_prime(Lk_int)
    if not prime:
        continue
    mod3 = Lk_int % 3
    eis = len(find_all_eis_norms(Lk_int, limit=200)) > 0
    gauss = len(find_all_gauss_norms(Lk_int)) > 0
    mass = m_e * PHI**k

    sm_match = ""
    sm_masses = {'e':0.000511,'mu':0.10566,'tau':1.77686,
                 'u':0.00216,'d':0.00467,'s':0.0934,
                 'c':1.275,'b':4.18,'t':172.76}
    for fname, fmass in sm_masses.items():
        err = abs(mass-fmass)/fmass*100
        if err < 8:
            sm_match = f"{fname}({err:.1f}%)"

    e_mark = "✓" if eis else " "
    g_mark = "✓" if gauss else " "
    print(f"{k:>4} {Lk_int:>10} {str(prime):>6} {mod3:>5} "
          f"{e_mark:>8} {g_mark:>10} {mass:>10.4f} {sm_match:>15}")
    all_results.append((k, Lk_int, prime, mod3, eis, gauss, mass, sm_match))

# ── CKM trace field via cusped manifold ──────────────────────────────────────
print()
print("="*65)
print("CKM TRACE FIELD: cusped ancestor approach")
print("="*65)
print()
print("The closed CKM manifold = m006(-5,2) Dehn filling of cusped m006")
print("The cusped m006 trace field IS computable via SnapPy")
print()

try:
    # Use the cusped manifold
    M_cusp_PMNS = snappy.Manifold("m003")
    M_cusp_CKM  = snappy.Manifold("m006")

    for label, M in [("PMNS cusped (m003)", M_cusp_PMNS),
                     ("CKM  cusped (m006)", M_cusp_CKM)]:
        print(f"{label}:")
        print(f"  vol = {float(M.volume()):.6f}")
        try:
            tf = M.invariant_trace_field()
            print(f"  Trace field: {tf}")
        except Exception as e:
            print(f"  Trace field error: {e}")
        try:
            qa = M.invariant_quaternion_algebra()
            print(f"  Quaternion algebra: {qa}")
        except Exception as e:
            print(f"  Quaternion algebra error: {e}")
        try:
            tf2 = M.trace_field()
            print(f"  Trace field (alt): {tf2}")
        except Exception as e:
            print(f"  Trace field (alt) error: {e}")
        print()
except Exception as e:
    print(f"Error: {e}")

# ── Gaussian norm Lucas primes (for Q(sqrt(-1)) case) ──────────────────────
print("="*65)
print("GAUSSIAN NORM LUCAS PRIMES (if CKM ~ Q(sqrt(-1)))")
print("="*65)
print()
print("A prime p is a Gaussian norm iff p ≡ 1 (mod 4) or p=2")
print()
print(f"{'k':>4} {'L_k_int':>10} {'mod4':>6} {'Gaussian?':>10} "
      f"{'witness':>12} {'mass_GeV':>10} {'SM quark?':>12}")
print("-"*65)

sm_quarks = {'u':0.00216,'d':0.00467,'s':0.0934,
             'c':1.275,'b':4.18,'t':172.76}

for k in range(0, 40):
    Lk = PHI**k + PHI**(-k)
    Lk_int = round(Lk)
    if abs(Lk - Lk_int) > 0.015: continue
    if Lk_int < 2: continue
    if not is_prime(Lk_int): continue

    mod4 = Lk_int % 4
    gauss_list = find_all_gauss_norms(Lk_int)
    is_gauss = len(gauss_list) > 0
    mass = m_e * PHI**k

    sm_match = ""
    for fname, fmass in sm_quarks.items():
        err = abs(mass-fmass)/fmass*100
        if err < 10:
            sm_match = f"{fname}({err:.1f}%)"

    if is_gauss or sm_match:
        w_str = str(gauss_list[0]) if gauss_list else "---"
        print(f"{k:>4} {Lk_int:>10} {mod4:>6} {str(is_gauss):>10} "
              f"{w_str:>12} {mass:>10.4f} {sm_match:>12}")

print()
print("="*65)
print("THE CENTRAL QUESTION")
print("="*65)
print()
print("PMNS sector: lepton masses ~ Lucas primes with Eisenstein norms")
print("  Trace field Q(sqrt(-3)), norm a^2 - ab + b^2")
print("  Eisenstein primes: p ≡ 1 (mod 3)")
print("  Confirmed: L_11=199 (mu), L_17=3571 (tau)")
print()
print("CKM sector: quark masses ~ Lucas primes with ??? norms")
print("  Trace field Q(??), norm ???")
print("  The quark masses have larger phi-lattice residuals (worse fit)")
print("  This is consistent with a different norm form")
print()
print("Key arithmetic fact:")
print("  199 ≡ 1 (mod 3) AND 199 ≡ 3 (mod 4) → Eisenstein YES, Gaussian NO")
print("  3571 ≡ 1 (mod 3) AND 3571 ≡ 3 (mod 4) → Eisenstein YES, Gaussian NO")
print()
print(f"  199 mod 4 = {199%4}")
print(f"  3571 mod 4 = {3571%4}")
print()
print("The lepton masses correspond to EISENSTEIN primes that are")
print("NOT Gaussian primes. This is the arithmetic signature of")
print("Q(sqrt(-3)) vs Q(sqrt(-1)).")
print()
print("If the quark masses correspond to GAUSSIAN primes, they would")
print("need to be ≡ 1 (mod 4). Let's check the quark mass ratios:")
for fname, fmass in sm_quarks.items():
    ratio = fmass/m_e
    # Find nearest Lucas prime
    best = min([(abs(PHI**k+PHI**(-k)-ratio), k, round(PHI**k+PHI**(-k)))
                for k in range(120)], key=lambda x: x[0])
    dist, k, Lk = best
    err = dist/ratio*100
    mod3 = Lk%3 if Lk>0 else -1
    mod4 = Lk%4 if Lk>0 else -1
    prime = is_prime(Lk) if 0 < Lk < 10**7 else False
    print(f"  {fname}: m/m_e={ratio:.1f}  "
          f"nearest_L_k={Lk}(k={k})  err={err:.2f}%  "
          f"mod3={mod3}  mod4={mod4}  prime={prime}")
