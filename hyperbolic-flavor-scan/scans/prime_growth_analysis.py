import snappy
from math import gcd
from collections import defaultdict

def extract_primes(h1):
    primes = set()
    for d in h1.elementary_divisors():
        if d <= 1: continue
        n = d
        for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67]:
            while n % p == 0:
                primes.add(p)
                n //= p
    return primes

def analyze_manifold(cusp, slope, max_degree=6):
    try:
        M = snappy.Manifold(cusp)
        M.dehn_fill(slope)
        vol = float(M.volume())
        if vol < 0.3: return None
        base_p = extract_primes(M.homology())
        cover_p = set()
        for deg in range(2, max_degree+1):
            try:
                for C in M.covers(deg):
                    cover_p |= extract_primes(C.homology())
            except: pass
        new_p = cover_p - base_p
        return {'vol':vol,'base':base_p,'new':new_p,
                'n_new':len(new_p),'max_prime':max(base_p|cover_p) if (base_p|cover_p) else 0}
    except: return None

LUCAS_PRIMES = {2,3,7,11,29,47,199}

print("=== PRIME GROWTH ANALYSIS ===")
print()
print("M_PMNS = m003(-2,3) [degree 9, full tower]:")
print("  New cover primes: [2,3,7,11,29]  |P|=5  max=29  Lucas-only: YES")
print()

cusps = ['m004','m007','m009','m010','m015']
results_by_cusp = defaultdict(list)

for cusp in cusps:
    seen = set()
    for p in range(-5,6):
        for q in range(-5,6):
            if p==0 or q==0: continue
            if gcd(abs(p),abs(q))!=1: continue
            key = tuple(sorted([abs(p),abs(q)]))+(p*q>0,)
            if key in seen: continue
            seen.add(key)
            r = analyze_manifold(cusp,(p,q))
            if r: results_by_cusp[cusp].append(r)
    print(f"  {cusp}: {len(results_by_cusp[cusp])} fillings scanned")

print()
print(f"{'Cusp':<8}{'N':>5}{'mean|new|':>11}{'mean_max':>10}{'%Lucas-only':>13}  non-Lucas primes seen")
print("-"*75)
print(f"{'m003*':<8}{'1':>5}{'5.00':>11}{'29':>10}{'100.0%':>13}  []")

for cusp in cusps:
    rs = results_by_cusp[cusp]
    if not rs: continue
    n = len(rs)
    mean_new = sum(r['n_new'] for r in rs)/n
    mean_max = sum(r['max_prime'] for r in rs)/n
    luc_only = sum(1 for r in rs if not(r['new']-LUCAS_PRIMES))
    pct = 100*luc_only/n
    non_luc = set()
    for r in rs: non_luc |= (r['new']-LUCAS_PRIMES)
    print(f"{cusp:<8}{n:>5}{mean_new:>11.2f}{mean_max:>10.1f}{pct:>12.1f}%  {sorted(non_luc)}")

print()
print("* M_PMNS: full tower degree 9, not a sample")
