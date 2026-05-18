# hecke_lucas_v2.sage
# Run: sage hecke_lucas_v2.sage
from sage.all import *

PHI = (1 + sqrt(5))/2
LOG_PHI = float(log(PHI))

def lucas_float(k):
    return float(PHI**k + PHI**(-k))

L = {}
for k in range(21):
    try:
        L[k] = ZZ(expand(PHI**k + PHI**(-k)))
    except:
        L[k] = lucas_float(k)

print("=" * 65)
print("Lucas numbers L_k = phi^k + phi^{-k}")
print("=" * 65)
for k in range(12):
    print(f"  L_{k:2d} = {L[k]}")
print()

print("Quarter-Lucas mass lattice positions:")
for q in [12, 18, 43, 44, 65, 68, 75, 106]:
    ell = q/4 * LOG_PHI
    print(f"  q={q:4d}: ell={ell:.5f}  phi^(q/4)={float(PHI**(q/4)):.4f}")
print()

print("=" * 65)
print("HECKE EIGENVALUES OF MODULAR FORMS")
print("=" * 65)
print()

for level, weight in [(3,2),(3,4),(5,2),(6,2),(15,2)]:
    try:
        S = Newforms(level, weight=weight, names='a')
        if not S:
            continue
        print(f"Level {level}, weight {weight}: {len(S)} newform(s)")
        for f in S[:2]:
            print(f"  Form: {f.label()}")
            print(f"  {'p':>4}  {'a_p':>10}  {'|a_p|':>8}  "
                  f"{'log|a_p|/logphi':>16}  {'nearest L_k':>14}")
            for p in [2,3,5,7,11,13,17,19,23,29]:
                try:
                    ap = f.hecke_eigenvalue(p)
                    ap_f = float(ap)
                    ratio = log(abs(ap_f))/LOG_PHI if abs(ap_f) > 0.001 else 0
                    nearest_k = min(range(20),
                                   key=lambda k: abs(lucas_float(k)-abs(ap_f)))
                    nearest_L = lucas_float(nearest_k)
                    dist = abs(abs(ap_f) - nearest_L)
                    match = "LUCAS!" if dist < 0.5 else ""
                    print(f"  {p:>4}  {ap_f:>10.4f}  {abs(ap_f):>8.4f}  "
                          f"{float(ratio):>16.6f}  "
                          f"L_{nearest_k}={nearest_L:.3f} {match}")
                except Exception as e:
                    print(f"  {p:>4}  error: {e}")
        print()
    except Exception as e:
        print(f"Level {level}, weight {weight}: {e}")
        print()

print("=" * 65)
print("LEVEL 5 NEWFORM (H1=Z/5 connection)")
print("=" * 65)
print()
try:
    S5 = Newforms(5, weight=2, names='a')
    if S5:
        f5 = S5[0]
        print(f"Newform: {f5.label()}")
        print(f"q-expansion: {f5.q_expansion(15)}")
        print()
        print(f"{'p':>4}  {'a_p':>8}  {'L_k match?':>12}")
        print("-"*35)
        for p in [2,3,5,7,11,13,17,19,23,29]:
            try:
                ap = ZZ(f5.hecke_eigenvalue(p))
                matches = [(k, L[k]) for k in range(15)
                           if abs(ap - L[k]) <= 1]
                mstr = f"≈ L_{matches[0][0]}={matches[0][1]}" if matches else ""
                print(f"{p:>4}  {ap:>8}  {mstr}")
            except Exception as e:
                print(f"{p:>4}  error: {e}")
except Exception as e:
    print(f"Error: {e}")

print()
print("=" * 65)
print("BRANDT MODULE (Quaternion algebra)")
print("=" * 65)
print()

for disc in [2, 3, 6, 10, 15]:
    try:
        A = QuaternionAlgebra(disc)
        O = A.maximal_order()
        B = BrandtModule(A, 1)
        print(f"Discriminant {disc}: rank={B.rank()}, "
              f"ramified at {A.ramified_primes()}")
        for p in [2, 3, 5, 7, 11, 13, 29]:
            try:
                T = B.hecke_matrix(p)
                evals = [float(e) for e in T.eigenvalues()]
                lucas_hits = []
                for e in evals:
                    for k in range(15):
                        if abs(abs(e) - lucas_float(k)) < 0.5:
                            lucas_hits.append(f"|{e:.3f}|≈L_{k}")
                hit_str = ", ".join(lucas_hits) if lucas_hits else ""
                print(f"  T_{p:2d} evals: {[f'{e:.3f}' for e in evals[:4]]}  "
                      f"{hit_str}")
            except Exception as e2:
                print(f"  T_{p}: {e2}")
        print()
        break  # use first working discriminant
    except Exception as e:
        print(f"Discriminant {disc}: {e}")

print()
print("=" * 65)
print("CONJECTURE SUMMARY")
print("=" * 65)
print()
print("Hecke-mass conjecture: a_p = L_{k(p)} * sqrt(N(p))")
print("where k(p) = degree at which p first appears in covering tower")
print()
covering = {2:9, 3:9, 7:9, 11:2, 29:5}
print(f"{'p':>4}  {'k':>4}  {'L_k':>8}  {'sqrt(N(p))':>12}  {'pred a_p':>10}")
print("-"*45)
for p, deg in covering.items():
    Lk = lucas_float(deg)
    Np = p if (p == 3 or p % 3 == 1) else p**2
    pred = Lk * Np**0.5
    print(f"{p:>4}  {deg:>4}  {Lk:>8.4f}  {Np**0.5:>12.4f}  {pred:>10.4f}")
print()
print("Compare predicted a_p values with Brandt/newform eigenvalues above.")
