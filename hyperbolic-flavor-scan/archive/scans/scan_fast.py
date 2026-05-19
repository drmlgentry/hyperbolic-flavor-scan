import snappy
import numpy as np
from math import gcd
import time

phi = (1+np.sqrt(5))/2
log_phi = np.log(phi)

def get_real_lengths(M, cutoff=4.0):
    spec = M.length_spectrum(cutoff)
    return [float(item.length.real()) for item in spec]

def log_n_hits(M, ns, cutoff=4.0, tol=0.02):
    try:
        lengths = get_real_lengths(M, cutoff)
        return [n for n in ns if any(abs(l - np.log(n)) < tol for l in lengths)]
    except:
        return []

def scan_cusp(cusp_name, p_range, q_range, max_abs=5):
    slopes = []
    for p in range(-max_abs, max_abs+1):
        for q in range(-max_abs, max_abs+1):
            if p == 0 or q == 0: continue
            if gcd(abs(p), abs(q)) != 1: continue
            # deduplicate (p,q) and (-p,-q)
            if (p,q) in slopes or (-p,-q) in slopes: continue
            slopes.append((p,q))
    
    ns = list(range(2, 61))
    results = {}
    total = len(slopes)
    
    print(f"\n=== Scanning {cusp_name} ({total} slopes) ===")
    start_time = time.time()
    
    for idx, (p,q) in enumerate(slopes):
        try:
            M = snappy.Manifold(cusp_name)
            M.dehn_fill((p,q))
            vol = float(M.volume())
            if vol < 0.5:
                print(f"  [{idx+1}/{total}] ({p},{q}): non-hyperbolic (vol={vol:.3f})")
                continue
            hits = log_n_hits(M, ns, cutoff=4.0, tol=0.02)
            results[(p,q)] = hits
            print(f"  [{idx+1}/{total}] ({p},{q}): vol={vol:.4f}, {len(hits)} hits")
        except Exception as e:
            print(f"  [{idx+1}/{total}] ({p},{q}): error: {e}")
        
        # progress every 10 items
        if (idx+1) % 10 == 0:
            elapsed = time.time() - start_time
            rate = (idx+1) / elapsed
            eta = (total - (idx+1)) / rate
            print(f"    -> progress: {idx+1}/{total} ({100*(idx+1)/total:.1f}%), ETA {eta:.0f}s")
    
    return results, slopes

# Run scans
print("="*60)
print("FAST SCAN: Dehn fillings of m003 and m006")
print("cutoff=4.0, tolerance=0.02, deduplicated slopes")
print("="*60)

m003_results, m003_slopes = scan_cusp('m003', range(-5,6), range(-5,6), max_abs=5)
m006_results, m006_slopes = scan_cusp('m006', range(-5,6), range(-5,6), max_abs=5)

print("\n" + "="*60)
print("SCAN COMPLETE")
print(f"m003: {len(m003_results)} hyperbolic fillings")
print(f"m006: {len(m006_results)} hyperbolic fillings")
print("="*60)

# Compute occurrence frequencies
ns = list(range(2, 61))
occurrence_m003 = {n: 0 for n in ns}
occurrence_m006 = {n: 0 for n in ns}

for hits in m003_results.values():
    for n in hits:
        occurrence_m003[n] += 1
for hits in m006_results.values():
    for n in hits:
        occurrence_m006[n] += 1

N_m003 = len(m003_results)
N_m006 = len(m006_results)

print("\n=== m003 OCCURRENCE FREQUENCY ===")
print(f"  Based on {N_m003} fillings")
print("  n   count     %")
for n in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]:
    cnt = occurrence_m003[n]
    pct = 100*cnt/N_m003
    star = "***" if cnt == N_m003 else "**" if cnt > 0.8*N_m003 else "*" if cnt > 0.5*N_m003 else ""
    print(f"  {n:3d}   {cnt:2d}/{N_m003:<2}  {pct:5.1f}%  {star}")

print("\n=== m006 OCCURRENCE FREQUENCY ===")
print(f"  Based on {N_m006} fillings")
print("  n   count     %")
for n in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]:
    cnt = occurrence_m006[n]
    pct = 100*cnt/N_m006
    star = "***" if cnt == N_m006 else "**" if cnt > 0.8*N_m006 else "*" if cnt > 0.5*N_m006 else ""
    print(f"  {n:3d}   {cnt:2d}/{N_m006:<2}  {pct:5.1f}%  {star}")

