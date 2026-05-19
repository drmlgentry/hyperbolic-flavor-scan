import snappy
from math import gcd
from collections import defaultdict

def extract_primes(h1):
    primes = set()
    for d in h1.elementary_divisors():
        if d <= 1: continue
        n = d
        for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]:
            while n % p == 0:
                primes.add(p)
                n //= p
    return primes

LUCAS_PRIMES = {2,3,7,11,29,47,199}

print("=== M_CKM = m006(-5,2) COVERING TOWER ===")
M_ckm = snappy.OrientableClosedCensus[43]
print(f"vol={float(M_ckm.volume()):.5f}  H1={M_ckm.homology()}")
base_p = extract_primes(M_ckm.homology())
print(f"Base primes: {sorted(base_p)}")

all_new = set()
for deg in range(2, 8):
    try:
        covers = M_ckm.covers(deg)
        deg_p = set()
        for C in covers: deg_p |= extract_primes(C.homology())
        new = deg_p - base_p - all_new
        all_new |= (deg_p - base_p)
        nl = new - LUCAS_PRIMES
        print(f"  deg {deg}: {len(covers)} covers  new={sorted(new)}  "
              f"non-Lucas={sorted(nl) if nl else '[]'}")
    except Exception as e:
        print(f"  deg {deg}: {e}")

print(f"\nAll new primes: {sorted(all_new)}")
print(f"Non-Lucas:      {sorted(all_new - LUCAS_PRIMES)}")
print(f"Lucas-pure:     {not bool(all_new - LUCAS_PRIMES)}")

print("\n=== m006 FILLING SCAN ===")
results = []
seen = set()
for p in range(-5,6):
    for q in range(-5,6):
        if p==0 or q==0: continue
        if gcd(abs(p),abs(q))!=1: continue
        key = tuple(sorted([abs(p),abs(q)]))+(p*q>0,)
        if key in seen: continue
        seen.add(key)
        try:
            M = snappy.Manifold('m006')
            M.dehn_fill((p,q))
            vol = float(M.volume())
            if vol < 0.3: continue
            bp = extract_primes(M.homology())
            cp = set()
            for deg in range(2,7):
                try:
                    for C in M.covers(deg): cp |= extract_primes(C.homology())
                except: pass
            new_p = cp - bp
            results.append({'vol':vol,'new':new_p,'n_new':len(new_p),
                           'max_p':max(bp|cp) if (bp|cp) else 0,
                           'lucas_pure':not bool(new_p-LUCAS_PRIMES)})
        except: pass

n = len(results)
luc = sum(1 for r in results if r['lucas_pure'])
mean_new = sum(r['n_new'] for r in results)/n
mean_max = sum(r['max_p'] for r in results)/n
all_nl = set()
for r in results: all_nl |= (r['new']-LUCAS_PRIMES)

print(f"m006: {n} fillings  mean|new|={mean_new:.2f}  "
      f"mean_max={mean_max:.1f}  %Lucas={100*luc/n:.1f}%")
print(f"Non-Lucas primes: {sorted(all_nl)}")
