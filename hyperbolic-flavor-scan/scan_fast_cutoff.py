import snappy
import numpy as np
from math import gcd
import time
import sys

phi = (1 + np.sqrt(5)) / 2
log_phi = np.log(phi)

BLACKLIST = [(-5, 4), (5, -4), (-4, 5), (4, -5)]

def is_blacklisted(p, q):
    return (p, q) in BLACKLIST or (-p, -q) in BLACKLIST

def scan_slope(p, q, cutoff=4.0, tol=0.02):
    try:
        M = snappy.Manifold('m003')
        M.dehn_fill((p, q))
        vol = float(M.volume())
        if vol < 0.5:
            return None, vol, 0
        spec = M.length_spectrum(cutoff)
        lengths = [float(x.length.real()) for x in spec]
        ns = list(range(2, 55))
        hits = [n for n in ns if any(abs(l - np.log(n)) < tol for l in lengths)]
        return hits, vol, len(lengths)
    except Exception as e:
        return None, -1, 0

def main():
    slopes = []
    seen = set()
    for p in range(-5, 6):
        for q in range(-5, 6):
            if p == 0 or q == 0: continue
            if gcd(abs(p), abs(q)) != 1: continue
            key = min((p, q), (-p, -q))
            if key in seen: continue
            seen.add(key)
            slopes.append((p, q))
    
    total = len(slopes)
    print(f"Total slopes: {total}")
    print(f"Blacklisted: {BLACKLIST}")
    print(f"cutoff=4.0, tol=0.02")
    print("=" * 60)
    sys.stdout.flush()
    
    hits_total = {n: 0 for n in range(2, 55)}
    hyp_count = 0
    
    for idx, (p, q) in enumerate(slopes):
        if is_blacklisted(p, q):
            print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) SKIPPED")
            sys.stdout.flush()
            continue
        
        start = time.time()
        hits, vol, n_len = scan_slope(p, q, cutoff=4.0, tol=0.02)
        elapsed = time.time() - start
        
        if hits is None:
            if vol < 0.5:
                print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) non-hyperbolic vol={vol:.3f} ({elapsed:.1f}s)")
            else:
                print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) ERROR ({elapsed:.1f}s)")
        else:
            hyp_count += 1
            for n in hits:
                hits_total[n] += 1
            print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) vol={vol:.4f} hits={len(hits):2d} geod={n_len:3d} ({elapsed:.1f}s)  hyp={hyp_count}")
        sys.stdout.flush()
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print(f"Hyperbolic fillings: {hyp_count}")
    print("\nOccurrence frequencies:")
    for n in [2, 3, 5, 7, 11, 13, 17, 18, 19, 21, 22, 23, 24, 25, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]:
        cnt = hits_total.get(n, 0)
        pct = 100 * cnt / hyp_count if hyp_count > 0 else 0
        print(f"  n={n:2d}: {cnt:2d}/{hyp_count} ({pct:5.1f}%)")

if __name__ == "__main__":
    main()
