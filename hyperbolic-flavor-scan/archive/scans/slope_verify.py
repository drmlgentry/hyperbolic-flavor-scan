import snappy, math
from math import gcd

PHI = (1+math.sqrt(5))/2
LOG_PHI = math.log(PHI)

print("=== SLOPE ENCODING VERIFICATION ===")
print()

print("Computing all m003 filling spectra...")
all_fillings = []
for p in range(-6,7):
    for q in range(-6,7):
        if p==0 or q==0: continue
        if gcd(abs(p),abs(q)) != 1: continue
        M = snappy.Manifold('m003')
        M.dehn_fill((p,q))
        try:
            vol = float(M.volume())
            if vol < 0.5: continue
            spec = M.length_spectrum(cutoff=3.5)
            geo = sorted([float(g.length.real()) for g in spec])
            all_fillings.append(((p,q), vol, geo))
        except: pass

print(f"Got {len(all_fillings)} fillings")
print()

def hits_log_n(geo, n, threshold=2.0):
    target = math.log(n)
    nearest = min(geo, key=lambda g: abs(g-target))
    delta = 100*abs(nearest-target)/target
    return delta < threshold, delta

print("Testing: do fillings m003(p,q) with |p|=k hit log(k) more than others?")
print()

for target_p in [2, 3, 4, 5, 6, 7]:
    with_p    = [(s,v,g) for s,v,g in all_fillings if abs(s[0])==target_p]
    without_p = [(s,v,g) for s,v,g in all_fillings if abs(s[0])!=target_p]

    hw  = [(s, *hits_log_n(g, target_p)) for s,v,g in with_p]
    hwo = [(s, *hits_log_n(g, target_p)) for s,v,g in without_p]

    n_w  = len(hw);  h_w  = sum(1 for _,hit,_ in hw  if hit)
    n_wo = len(hwo); h_wo = sum(1 for _,hit,_ in hwo if hit)
    pct_w  = 100*h_w/n_w   if n_w  > 0 else 0
    pct_wo = 100*h_wo/n_wo if n_wo > 0 else 0

    best_w  = sorted([(d,s) for s,h,d in hw],  key=lambda x:x[0])[:5]
    best_wo = sorted([(d,s) for s,h,d in hwo], key=lambda x:x[0])[:3]

    ratio_str = f"{pct_w/pct_wo:.2f}x" if pct_wo > 0 else "inf (0% without)"
    print(f"|p|={target_p}: log({target_p})={math.log(target_p):.4f}")
    print(f"  WITH    |p|={target_p}: {h_w:>3}/{n_w:<3} = {pct_w:>5.1f}%")
    print(f"  WITHOUT |p|={target_p}: {h_wo:>3}/{n_wo:<3} = {pct_wo:>5.1f}%")
    print(f"  Enrichment ratio: {ratio_str}")
    print(f"  Best WITH hits:")
    for delta, slope in best_w[:5]:
        print(f"    m003{slope}: {delta:.4f}%")
    print(f"  Best WITHOUT hits:")
    for delta, slope in best_wo[:3]:
        print(f"    m003{slope}: {delta:.4f}%")
    print()
