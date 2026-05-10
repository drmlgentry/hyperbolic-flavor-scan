import snappy
import numpy as np
from math import gcd
import time
import sys
import os

phi = (1 + np.sqrt(5)) / 2
log_phi = np.log(phi)

def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')

def get_real_lengths(M, cutoff=5.0):
    spec = M.length_spectrum(cutoff)
    return [float(item.length.real()) for item in spec]

def log_n_hits(M, ns, cutoff=5.0, tol=0.02):
    try:
        lengths = get_real_lengths(M, cutoff)
        return [n for n in ns if any(abs(l - np.log(n)) < tol for l in lengths)]
    except:
        return []

def print_table(occ, total, label):
    print(f"\n{label} ({total} fillings processed so far)")
    print("   n   count     %  note")
    for n in [2,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]:
        cnt = occ.get(n, 0)
        pct = 100 * cnt / total if total > 0 else 0
        star = "***" if cnt == total else "**" if cnt > 0.8*total else "*" if cnt > 0.5*total else ""
        print(f"  {n:3d}   {cnt:2d}/{total:<2}  {pct:5.1f}%  {star}")

def scan_cusp(cusp_name, max_abs=5, cutoff=5.0, tol=0.02):
    slopes = []
    for p in range(-max_abs, max_abs+1):
        for q in range(-max_abs, max_abs+1):
            if p == 0 or q == 0: continue
            if gcd(abs(p), abs(q)) != 1: continue
            if (p,q) in slopes or (-p,-q) in slopes: continue
            slopes.append((p,q))
    
    ns = list(range(2, 61))
    results = {}
    occ = {n: 0 for n in ns}
    hyperbolic = []
    total = len(slopes)
    
    print(f"\n{'='*60}")
    print(f"LIVE SCAN: {cusp_name} (cutoff={cutoff}, tol={tol})")
    print(f"Total slopes: {total}")
    print(f"{'='*60}")
    
    start_time = time.time()
    for idx, (p,q) in enumerate(slopes):
        try:
            M = snappy.Manifold(cusp_name)
            M.dehn_fill((p,q))
            vol = float(M.volume())
            if vol < 0.5:
                continue
            hits = log_n_hits(M, ns, cutoff=cutoff, tol=tol)
            results[(p,q)] = hits
            hyperbolic.append((p,q,vol,len(hits)))
            for n in hits:
                occ[n] += 1
            
            # Update display every 5 fillings
            if (idx+1) % 5 == 0 or idx+1 == total:
                clear_screen()
                print(f"\n{cusp_name}: {len(hyperbolic)} hyperbolic so far / {total} slopes")
                print(f"Progress: {idx+1}/{total} ({100*(idx+1)/total:.1f}%)")
                elapsed = time.time() - start_time
                rate = (idx+1) / elapsed if elapsed > 0 else 0
                eta = (total - (idx+1)) / rate if rate > 0 else 0
                print(f"ETA: {eta:.0f}s")
                print_table(occ, len(hyperbolic), cusp_name.upper())
                
        except Exception as e:
            continue
    
    clear_screen()
    print(f"\n{'='*60}")
    print(f"{cusp_name} COMPLETE")
    print(f"Hyperbolic fillings: {len(hyperbolic)} out of {total}")
    print(f"{'='*60}")
    print_table(occ, len(hyperbolic), cusp_name.upper())
    
    return results, hyperbolic, occ

# Run scans
m003_results, m003_hyp, m003_occ = scan_cusp('m003', max_abs=5, cutoff=5.0, tol=0.02)
m006_results, m006_hyp, m006_occ = scan_cusp('m006', max_abs=5, cutoff=5.0, tol=0.02)

# Final summary
print("\n\n" + "="*60)
print("FINAL SUMMARY")
print("="*60)
print(f"m003: {len(m003_hyp)} hyperbolic fillings")
print(f"m006: {len(m006_hyp)} hyperbolic fillings")
print("\nThreshold n (>80% occurrence):")
for label, occ, N in [("m003", m003_occ, len(m003_hyp)), ("m006", m006_occ, len(m006_hyp))]:
    thr = min([n for n in range(2,61) if occ.get(n,0) > 0.8*N], default=None)
    print(f"  {label}: n={thr}  (log(n)/log(phi)={np.log(thr)/log_phi:.4f})")

