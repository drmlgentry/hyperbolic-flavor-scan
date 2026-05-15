"""
Trace field analysis and Lucas-Eisenstein norm structure
for M_PMNS = m003(-2,3) and M_CKM = m006(-5,2)

Key question: does the Lucas-Eisenstein theorem extend to the
quark sector via the trace field of M_CKM?
"""
import snappy, numpy as np
import warnings; warnings.filterwarnings('ignore')

PHI = (1+5**0.5)/2
LOG_PHI = np.log(PHI)
m_e = 0.000511  # GeV

def is_prime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5)+1, 2):
        if n%i==0: return False
    return True

def eis_norm(a, b):
    """Norm in Z[omega]: N(a+b*omega) = a^2 - ab + b^2"""
    return a*a - a*b + b*b

def is_norm_in_Zomega(N, limit=200):
    """Check if N is representable as a^2 - ab + b^2"""
    for a in range(-limit, limit+1):
        disc = 4*N - 3*a*a
        if disc < 0: continue
        sq = int(disc**0.5 + 0.5)
        for s in [sq-1, sq, sq+1]:
            if s >= 0 and s*s == disc:
                for sign in [1, -1]:
                    num = a + sign*s
                    if num % 2 == 0:
                        b = num // 2
                        if eis_norm(a,b) == N:
                            return True, (a,b)
    return False, None

print("="*70)
print("TRACE FIELD VERIFICATION: M_PMNS and M_CKM")
print("="*70)
print()

M_PMNS = snappy.OrientableClosedCensus[1]
M_CKM  = snappy.OrientableClosedCensus[43]

for label, M in [("PMNS", M_PMNS), ("CKM", M_CKM)]:
    print(f"M_{label} = {M.name()}")
    print(f"  vol = {float(M.volume()):.6f}")
    print(f"  H1  = {M.homology()}")
    try:
        tf = M.invariant_trace_field()
        print(f"  Trace field: {tf}")
    except Exception as e:
        print(f"  Trace field: (error: {e})")
    try:
        qa = M.invariant_quaternion_algebra()
        print(f"  Quaternion algebra: {qa}")
    except Exception as e:
        print(f"  Quaternion algebra: (error: {e})")
    try:
        cs = float(M.chern_simons())
        print(f"  Chern-Simons: {cs:.6f}")
    except: pass
    print()

# ── Lucas-Eisenstein analysis for both manifolds ────────────────────────────
print("="*70)
print("LUCAS-EISENSTEIN NORM STRUCTURE")
print("="*70)
print()

# SM fermion masses and their phi-lattice quantum numbers
fermions = {
    'e':   (0.000511,   0),
    'mu':  (0.10566,   44),
    'tau': (1.77686,   68),
    'u':   (0.00216,   12),
    'd':   (0.00467,   18),
    's':   (0.0934,    43),
    'c':   (1.275,     65),
    'b':   (4.18,      75),
    't':   (172.76,   106),
}

print("SM fermion masses and Lucas numbers:")
print(f"  {'f':>4} {'q':>5} {'m(GeV)':>10} {'phi^(q/4)':>12} "
      f"{'nearest_L_k':>12} {'k':>4} {'error%':>8}")
print("  "+"-"*65)
for name, (mass, q) in fermions.items():
    if q == 0: continue
    pred_mass = m_e * PHI**(q/4)
    # Find nearest Lucas number to mass ratio m/m_e
    ratio = mass/m_e
    best_k = min(range(120), key=lambda k: abs(PHI**k + PHI**(-k) - ratio))
    Lk = PHI**best_k + PHI**(-best_k)
    err = abs(ratio - Lk)/ratio*100
    # Also check near integer
    Lk_int = round(Lk)
    err_int = abs(ratio - Lk_int)/ratio*100
    print(f"  {name:>4} {q:>5} {mass:>10.5f} {pred_mass:>12.4f} "
          f"{Lk_int:>12} {best_k:>4} {err_int:>8.3f}%")

print()

# ── The key check: which Lucas numbers are Eisenstein norms in Z[omega] ──────
print("="*70)
print("LUCAS PRIMES AS EISENSTEIN NORMS IN Z[omega] (trace field of M_PMNS)")
print("="*70)
print()
print("Z[omega] has norm form N(a+b*omega) = a^2 - ab + b^2")
print("A prime p is an Eisenstein norm iff p ≡ 1 (mod 3)")
print()
print(f"{'k':>4} {'L_k (int)':>12} {'prime':>6} {'mod3':>5} "
      f"{'norm?':>7} {'witness':>14} {'mass_GeV':>12} {'SM?':>15}")
print("  "+"-"*75)

predictions = []
for k in range(0, 50):
    Lk_exact = PHI**k + PHI**(-k)
    Lk_int = round(Lk_exact)
    if abs(Lk_exact - Lk_int) > 0.01:
        continue
    if Lk_int < 2:
        continue
    
    prime = is_prime(Lk_int)
    mod3 = Lk_int % 3
    
    if not prime or mod3 != 1:
        continue
    
    is_norm, witness = is_norm_in_Zomega(Lk_int, limit=150)
    mass_GeV = m_e * PHI**k
    
    # Find SM particle
    sm_match = ""
    best_err = 999
    for fname, (fmass, fq) in fermions.items():
        err = abs(mass_GeV - fmass)/fmass*100
        if err < best_err:
            best_err = err
            if err < 10:
                sm_match = f"{fname} ({err:.1f}%)"
    
    mark = "<<<" if (is_norm and best_err < 5) else ("<<" if is_norm else "")
    w_str = str(witness) if witness else "---"
    print(f"  {k:>4} {Lk_int:>12} {str(prime):>6} {mod3:>5} "
          f"{str(is_norm):>7} {w_str:>14} {mass_GeV:>12.4f} "
          f"{sm_match:>15} {mark}")
    
    if is_norm:
        predictions.append((k, Lk_int, mass_GeV, sm_match, witness))

print()
print("Summary of Lucas prime Eisenstein norms:")
print(f"  k = {[p[0] for p in predictions]}")
print(f"  L_k = {[p[1] for p in predictions]}")
print()

# ── Now check what norm form the CKM trace field has ────────────────────────
print("="*70)
print("M_CKM TRACE FIELD NORM ANALYSIS")
print("="*70)
print()
print("M_CKM has trace field Q(sqrt(-1)) or Q(sqrt(-3))?")
print("Need to check: what is the norm form of the CKM trace field?")
print()
print("If trace field is Q(sqrt(-d)), norm form is N(a+b*sqrt(-d)) = a^2 + d*b^2")
print("or for Z[omega] type: a^2 - ab + b^2")
print()

# The CKM trace field discriminant
# m006 has trace field Q(sqrt(-1)) based on prior work
# Let's verify by checking which norms appear as geodesic lengths

rho_CKM = M_CKM.polished_holonomy()
rho_PMNS = M_PMNS.polished_holonomy()

def get_trace(rho, word):
    try:
        mat = np.array(rho(word), dtype=complex)
        mat = mat/np.sqrt(np.linalg.det(mat))
        return float(np.real(np.trace(mat)))
    except: return None

print("Holonomy traces for short words (real part):")
print("These should lie in the trace field ring of integers")
print()

for label, rho, M in [("PMNS", rho_PMNS, M_PMNS),
                       ("CKM",  rho_CKM,  M_CKM)]:
    print(f"M_{label}:")
    print(f"  {'word':>10}  {'Re(tr)':>10}  {'round':>6}  "
          f"{'as a^2+b^2':>12}  {'as a^2-ab+b^2':>14}")
    for word in ['a','b','aa','ab','aB','aaB','baa',
                 'aabb','abab','aaBB']:
        tr = get_trace(rho, word)
        if tr is None: continue
        tr_int = round(tr)
        # Check if |tr_int| is a norm in Z[i] (a^2+b^2)
        # and in Z[omega] (a^2-ab+b^2)
        absval = abs(tr_int)
        
        # Z[i] norm: N(a+bi) = a^2+b^2
        is_Zi = any(a*a+b*b == absval 
                    for a in range(int(absval**0.5)+2)
                    for b in range(int(absval**0.5)+2))
        # Z[omega] norm
        is_Zo, _ = is_norm_in_Zomega(absval, limit=20) if absval > 0 else (False,None)
        
        print(f"  {word:>10}  {tr:>10.4f}  {tr_int:>6}  "
              f"{'YES' if is_Zi else 'no':>12}  "
              f"{'YES' if is_Zo else 'no':>14}")
    print()

# ── CKM quark mass prediction ────────────────────────────────────────────────
print("="*70)
print("QUARK MASS LUCAS-EISENSTEIN CHECK")
print("="*70)
print()
print("If M_CKM has trace field Q(sqrt(-1)), the relevant norm is")
print("the Gaussian integer norm: N(a+bi) = a^2 + b^2")
print()

def is_gaussian_norm(N, limit=200):
    """Check if N = a^2 + b^2"""
    for a in range(0, int(N**0.5)+1):
        b_sq = N - a*a
        b = int(b_sq**0.5+0.5)
        if b*b == b_sq:
            return True, (a,b)
    return False, None

print("Checking SM fermion mass ratios against Gaussian integer norms:")
print(f"  {'f':>4} {'ratio=m/m_e':>14} {'nearest L_k':>12} {'k':>4} "
      f"{'Gaussian norm?':>15} {'witness':>10}")
print("  "+"-"*65)
for name, (mass, q) in fermions.items():
    if q == 0: continue
    ratio = mass/m_e
    best_k = min(range(120), key=lambda k: abs(PHI**k+PHI**(-k)-ratio))
    Lk = round(PHI**best_k + PHI**(-best_k))
    err = abs(ratio-Lk)/ratio*100 if Lk > 0 else 999
    is_gn, witness_gn = is_gaussian_norm(Lk) if Lk > 0 and err < 10 else (False, None)
    if err < 10:
        print(f"  {name:>4} {ratio:>14.2f} {Lk:>12} {best_k:>4} "
              f"{str(is_gn):>15} {str(witness_gn):>10}")

print()
print("="*70)
print("FINAL: LUCAS-EISENSTEIN THEOREM SCOPE")
print("="*70)
print()
print("M_PMNS (trace field Q(sqrt(-3)), norm = a^2-ab+b^2):")
for k, Lk, mass, sm, witness in predictions:
    print(f"  k={k:2d}: L_k={Lk:<8} -> {mass:.4f} GeV  {sm}  witness={witness}")
print()
print("Open questions:")
print("  1. Does M_CKM (trace field Q(sqrt(-1))?) have analogous")
print("     Gaussian norm Lucas primes encoding quark masses?")
print("  2. Is k=19 (L_19=9349, 4.78 GeV) the b quark or BSM?")
print("  3. Can this be proven from the arithmetic of the quaternion")
print("     algebra over Q(sqrt(-3))?")
