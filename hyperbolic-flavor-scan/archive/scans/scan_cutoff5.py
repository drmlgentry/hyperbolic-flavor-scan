import snappy
import numpy as np
from math import gcd
import time
import sys

phi = (1 + np.sqrt(5)) / 2
log_phi = np.log(phi)

def get_real_lengths(M, cutoff=5.0):
    spec = M.length_spectrum(cutoff)
    return [float(item.length.real()) for item in spec]

def log_n_hits(M, ns, cutoff=5.0, tol=0.02):
    try:
        lengths = get_real_lengths(M, cutoff)
        return [n for n in ns if any(abs(l - np.log(n)) < tol for l in lengths)]
    except:
        return []

def progress_bar(current, total, start_time, label=""):
    elapsed = time.time() - start_time
    pct = 100 * current / total
    eta = (elapsed / current) * (total - current) if current > 0 else 0
    bar_len = 40
    filled = int(bar_len * current / total)
    bar = '█' * filled + '░' * (bar_len - filled)
    sys.stdout.write(f"\r{label} [{bar}] {current}/{total} ({pct:.1f}%) ETA {eta:.0f}s")
    sys.stdout.flush()

def scan_cusp(cusp_name, max_abs=5, cutoff=5.0, tol=0.02):
    slopes = []
    for p in range(-max_abs, max_abs+1):
        for q in range(-max_abs, max_abs+1):
            if p == 0 or q == 0:
                continue
            if gcd(abs(p), abs(q)) != 1:
                continue
            if (p, q) in slopes or (-p, -q) in slopes:
                continue
            slopes.append((p, q))
    
    ns = list(range(2, 61))
    results = {}
    hyperbolic = []
    total = len(slopes)
    
    print(f"\n{'='*60}")
    print(f"Scanning {cusp_name} (cutoff={cutoff}, tol={tol})")
    print(f"Slopes to test: {total}")
    print(f"{'='*60}")
    
    start_time = time.time()
    for idx, (p, q) in enumerate(slopes):
        try:
            M = snappy.Manifold(cusp_name)
            M.dehn_fill((p, q))
            vol = float(M.volume())
            if vol < 0.5:
                progress_bar(idx+1, total, start_time, f"{cusp_name} (skipping non-hyperbolic)")
                continue
            hits = log_n_hits(M, ns, cutoff=cutoff, tol=tol)
            results[(p, q)] = hits
            hyperbolic.append((p, q, vol, len(hits)))
            progress_bar(idx+1, total, start_time, f"{cusp_name}")
        except Exception as e:
            progress_bar(idx+1, total, start_time, f"{cusp_name} (error)")
            continue
    
    print(f"\n\n{cusp_name}: {len(hyperbolic)} hyperbolic fillings out of {total}")
    return results, hyperbolic

# Run scans
m003_results, m003_hyp = scan_cusp('m003', max_abs=5, cutoff=5.0, tol=0.02)
m006_results, m006_hyp = scan_cusp('m006', max_abs=5, cutoff=5.0, tol=0.02)

# Count occurrences
ns = list(range(2, 61))
occ_m003 = {n: 0 for n in ns}
occ_m006 = {n: 0 for n in ns}

for hits in m003_results.values():
    for n in hits:
        occ_m003[n] += 1
for hits in m006_results.values():
    for n in hits:
        occ_m006[n] += 1

N_m003 = len(m003_results)
N_m006 = len(m006_results)

print("\n" + "="*60)
print("SUMMARY")
print("="*60)

print(f"\nm003: {N_m003} hyperbolic fillings")
print("   n   count     %  note")
for n in [2,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]:
    cnt = occ_m003[n]
    pct = 100 * cnt / N_m003
    star = "***" if cnt == N_m003 else "**" if cnt > 0.8*N_m003 else "*" if cnt > 0.5*N_m003 else ""
    print(f"  {n:3d}   {cnt:2d}/{N_m003:<2}  {pct:5.1f}%  {star}")

print(f"\nm006: {N_m006} hyperbolic fillings")
print("   n   count     %  note")
for n in [2,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]:
    cnt = occ_m006[n]
    pct = 100 * cnt / N_m006
    star = "***" if cnt == N_m006 else "**" if cnt > 0.8*N_m006 else "*" if cnt > 0.5*N_m006 else ""
    print(f"  {n:3d}   {cnt:2d}/{N_m006:<2}  {pct:5.1f}%  {star}")

# Find threshold
threshold_m003 = min([n for n in ns if occ_m003[n] > 0.8 * N_m003], default=None)
threshold_m006 = min([n for n in ns if occ_m006[n] > 0.8 * N_m006], default=None)
print(f"\nThreshold n (80% occurrence): m003: {threshold_m003}, m006: {threshold_m006}")
print(f"log(threshold)/log(phi): m003: {np.log(threshold_m003)/log_phi:.4f}, m006: {np.log(threshold_m006)/log_phi:.4f}")

