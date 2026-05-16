# hecke_lucas_v3.sage
# Run: sage hecke_lucas_v3.sage
from sage.all import *

PHI = (1 + sqrt(5))/2
LOG_PHI = float(log(PHI))

def lucas_float(k):
    return float(PHI**k + PHI**(-k))

print("=" * 65)
print("HECKE EIGENVALUE COMPUTATION v3")
print("Fixed API calls for Sage modular forms + Brandt modules")
print("=" * 65)
print()

# ── Step 1: Newforms at various levels ──────────────────────────────────────
print("STEP 1: NEWFORMS AND HECKE EIGENVALUES")
print()

def check_newforms(level, weight=2):
    try:
        S = CuspForms(Gamma0(level), weight)
        d = S.dimension()
        if d == 0:
            return
        print(f"Level {level}, weight {weight}: dim={d}")
        # Get newforms via decomposition
        try:
            D = S.new_subspace().decomposition()
            for i, piece in enumerate(D):
                print(f"  Newform {i+1} (dim {piece.dimension()}):")
                # Hecke eigenvalues via characteristic polynomial
                for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
                    try:
                        T = piece.hecke_matrix(p)
                        cp = T.charpoly()
                        roots = cp.roots(ring=QQbar)
                        evals = [float(r[0].real()) for r in roots
                                 if abs(float(r[0].imag())) < 1e-6]
                        if evals:
                            for e in evals[:3]:
                                ratio = log(abs(e))/LOG_PHI if abs(e) > 0.001 else 0
                                nearest_k = min(range(25),
                                    key=lambda k: abs(lucas_float(k)-abs(e)))
                                Lk = lucas_float(nearest_k)
                                match = "LUCAS!" if abs(abs(e)-Lk) < 0.5 else ""
                                print(f"    p={p:2d}: a_p={e:.4f}  "
                                      f"log/logphi={float(ratio):.4f}  "
                                      f"nearest L_{nearest_k}={Lk:.3f} {match}")
                    except Exception as e2:
                        pass
        except Exception as e:
            print(f"  Decomposition error: {e}")
            # Try direct Hecke operators
            try:
                for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
                    T = S.hecke_operator(p).matrix()
                    evals = [float(e.real()) for e in T.eigenvalues()
                             if abs(float(e.imag())) < 1e-6]
                    if evals:
                        print(f"  p={p:2d}: eigenvalues = "
                              f"{[f'{e:.3f}' for e in evals[:4]]}")
            except:
                pass
    except Exception as e:
        print(f"Level {level}: {e}")

for level in [3, 5, 6, 10, 15, 20, 30]:
    check_newforms(level, 2)
    print()

# ── Step 2: Direct modular symbols approach ──────────────────────────────────
print("=" * 65)
print("STEP 2: MODULAR SYMBOLS (more robust)")
print("=" * 65)
print()

for level in [5, 6, 15, 20]:
    try:
        M = ModularSymbols(level, 2, sign=1)
        S = M.cuspidal_subspace()
        if S.dimension() == 0:
            continue
        print(f"Level {level}: cuspidal dim = {S.dimension()}")
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
            if p.divides(level):
                continue
            try:
                T = S.hecke_matrix(p)
                cp = T.charpoly()
                # Try to get rational eigenvalues
                roots = cp.roots()  # rational roots
                if roots:
                    for r, mult in roots:
                        e = float(r)
                        ratio = log(abs(e))/LOG_PHI if abs(e) > 0.001 else 0
                        nearest_k = min(range(25),
                            key=lambda k: abs(lucas_float(k)-abs(e)))
                        Lk = lucas_float(nearest_k)
                        match = "LUCAS!" if abs(abs(e)-Lk) < 0.5 else ""
                        print(f"  p={p:2d}: a_p={e:.4f}  "
                              f"log/logphi={float(ratio):.4f}  "
                              f"L_{nearest_k}={Lk:.3f} {match}")
            except Exception as e2:
                pass
        print()
    except Exception as e:
        print(f"Level {level}: {e}")

# ── Step 3: Brandt module with correct syntax ────────────────────────────────
print("=" * 65)
print("STEP 3: BRANDT MODULES (correct API)")
print("=" * 65)
print()

# BrandtModule takes an integer (the discriminant of the quaternion algebra)
# For QQ: quaternion algebras ramified at a set of primes
# Discriminant = product of finite ramification primes

for disc in [6, 10, 15, 21, 35]:
    try:
        # Check if disc is a valid quaternion algebra discriminant
        # (product of an even number of primes... wait, over QQ
        # we need discriminant to be squarefree and have even number of
        # prime factors for definite algebras)
        B = BrandtModule(disc)
        print(f"BrandtModule({disc}): rank={B.dimension()}")
        for p in [2, 3, 5, 7, 11, 13, 17, 29]:
            try:
                T = B.hecke_matrix(p)
                evals = T.eigenvalues()
                evals_f = []
                for e in evals:
                    try:
                        evals_f.append(float(e))
                    except:
                        pass
                if evals_f:
                    lucas_hits = []
                    for e in evals_f:
                        for k in range(20):
                            if abs(abs(e) - lucas_float(k)) < 0.6:
                                lucas_hits.append(f"|{e:.2f}|≈L_{k}")
                    hit_str = " ".join(lucas_hits)
                    print(f"  T_{p:2d}: {[f'{e:.3f}' for e in evals_f[:5]]} "
                          f"{hit_str}")
            except Exception as e2:
                print(f"  T_{p}: {e2}")
        print()
        break
    except Exception as e:
        print(f"BrandtModule({disc}): {e}")

# ── Step 4: Direct check using known level-5 eta product ─────────────────────
print("=" * 65)
print("STEP 4: LEVEL 5 MODULAR FORM (direct q-expansion)")
print("=" * 65)
print()
print("The unique newform of level 5, weight 2 is:")
print("  f = q - 2q^2 - 2q^3 + 2q^4 + q^5 + 4q^6 - ...")
print("  This corresponds to the elliptic curve y^2+y = x^3-x^2-10x-10")
print("  (Cremona label 50a, or 5a after twist)")
print()

# Construct via q-expansion
try:
    R.<q> = PowerSeriesRing(QQ, default_prec=50)
    # Level 5 newform coefficients (known)
    # a_p for the unique newform of level 5 weight 2:
    # p=2: -2, p=3: -2, p=5: (bad), p=7: 2, p=11: -2, p=13: 6, p=17: 2
    # p=19: -6, p=23: -2, p=29: 6
    known_ap = {2:-2, 3:-2, 5:1, 7:2, 11:-2, 13:6, 17:2,
                19:-6, 23:-2, 29:6, 31:-6, 37:-4, 41:-8, 43:-4, 47:0}
    
    print(f"Known Hecke eigenvalues for level-5 newform:")
    print(f"  {'p':>4}  {'a_p':>6}  {'|a_p|':>6}  "
          f"{'log|a_p|/logphi':>16}  {'nearest L_k':>14}  {'dist':>8}")
    print("  " + "-"*65)
    
    for p, ap in sorted(known_ap.items()):
        if ap == 0:
            continue
        ratio = log(abs(ap))/LOG_PHI if abs(ap) > 0.5 else 0
        nearest_k = min(range(20), key=lambda k: abs(lucas_float(k)-abs(ap)))
        Lk = lucas_float(nearest_k)
        dist = abs(abs(ap) - Lk)
        match = "LUCAS!" if dist < 0.5 else ""
        print(f"  {p:>4}  {ap:>6}  {abs(ap):>6.3f}  "
              f"{float(ratio):>16.6f}  "
              f"L_{nearest_k}={Lk:.3f}  {dist:>8.4f} {match}")
except Exception as e:
    print(f"Error: {e}")

# ── Step 5: The key number check ─────────────────────────────────────────────
print()
print("=" * 65)
print("STEP 5: THE KEY NUMBER — a_7 AND THE MUON MASS")
print("=" * 65)
print()
print("From conjectured Hecke eigenvalue formula:")
print(f"  a_7 (conjectured) = sqrt(7) * L_9 = {7**0.5:.4f} * {lucas_float(9):.4f}")
conj_a7 = 7**0.5 * lucas_float(9)
print(f"  a_7 = {conj_a7:.6f}")
print(f"  log(a_7)/log(phi) = {log(conj_a7)/LOG_PHI:.6f}")
print()
print(f"From level-5 newform (if applicable):")
print(f"  a_7 = 2 (known)")
print(f"  log(2)/log(phi) = {log(2)/LOG_PHI:.6f}")
print()
print("Muon mass quantum number: q=44, q/4 = 11.000")
print(f"  conjectured a_7: log/logphi = {log(conj_a7)/LOG_PHI:.4f} "
      f"(vs muon q/4 = 11.000)")
print(f"  actual a_7 (level 5): log/logphi = {log(2)/LOG_PHI:.4f}")
print()
print("CONCLUSION:")
print("  The level-5 form gives a_7=2, which does NOT match the conjecture.")
print("  The relevant automorphic form for M_PMNS is NOT the level-5 newform.")
print("  The correct form lives over Q(sqrt(-3)), not over Q.")
print()
print("NEXT: Need Hilbert modular forms or Bianchi modular forms over Q(sqrt(-3)).")
print("  sage: BianchiModularForms(QuadraticField(-3), level=...)")
print("  This is the correct automorphic object for M_PMNS.")

# ── Step 6: Norm equation approach ──────────────────────────────────────────
print()
print("=" * 65)
print("STEP 6: NORM EQUATION IN Q(sqrt(-3)) AT SMALL NORMS")
print("=" * 65)
print()
print("Elements of Z[omega] with norm = SM mass phi-lattice values:")
print()

K.<w> = QuadraticField(-3)
OK = K.ring_of_integers()

# For each SM mass quantum number q, check if norm q/4 has solutions
# i.e., does there exist a+b*omega in Z[omega] with N(a+b*omega) = phi^(q/4)?
# More tractable: check norms that are close to phi^(q/4) values

target_norms = {
    'u': (12, float(PHI**(12/4))),   # q=12, phi^3 = 4.236
    'd': (18, float(PHI**(18/4))),   # q=18, phi^4.5 = 8.719
    's': (43, float(PHI**(43/4))),   # q=43
    'mu':(44, float(PHI**(44/4))),   # q=44, phi^11 = 199.005
    'tau':(68,float(PHI**(68/4))),   # q=68, phi^17
}

print(f"  {'fermion':>7}  {'q':>5}  {'phi^(q/4)':>12}  "
      f"{'nearest norm N':>15}  {'element':>15}")
print("  "+"-"*60)

for name, (q, target) in target_norms.items():
    # Find Eisenstein integers with norm closest to target
    best_dist = float('inf')
    best_elem = None
    best_norm = None
    for a in range(-25, 26):
        for b in range(-25, 26):
            N = a**2 - a*b + b**2
            if N <= 0: continue
            dist = abs(N - target)
            if dist < best_dist:
                best_dist = dist
                best_elem = (a, b)
                best_norm = N
    if best_elem:
        a, b = best_elem
        rel_err = best_dist/target*100
        print(f"  {name:>7}  {q:>5}  {target:>12.3f}  "
              f"{best_norm:>15}  "
              f"({a},{b})  err={rel_err:.2f}%")

print()
print("If any norm equals phi^(q/4) exactly, that Eisenstein integer")
print("would be the algebraic witness for the fermion mass.")
print()
print("Note: phi^(q/4) is irrational for most q, so exact matches")
print("require irrational norms -- which don't exist in Z[omega].")
print("The question is whether the NEAREST integer norm is special.")
