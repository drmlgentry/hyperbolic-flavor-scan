import snappy
import numpy as np
from math import gcd
import time
import sys

phi = (1 + np.sqrt(5)) / 2
log_phi = np.log(phi)

# List of slopes known to be slow or problematic
SLOW_SLOPES = [
    (-5, 4), (5, -4),  # the one that hung
    # Add more as you encounter them:
    # (-5, -4), (5, 4),
    # (-4, 5), (4, -5),
]

def is_blacklisted(p, q):
    return (p, q) in SLOW_SLOPES or (-p, -q) in SLOW_SLOPES

def get_real_lengths(M, cutoff=4.5):
    try:
        spec = M.length_spectrum(cutoff)
        return [float(item.length.real()) for item in spec]
    except:
        return []

def log_n_hits(M, ns, cutoff=4.5, tol=0.02):
    try:
        lengths = get_real_lengths(M, cutoff)
        return [n for n in ns if any(abs(l - np.log(n)) < tol for l in lengths)]
    except:
        return []

def scan_cusp(cusp_name, max_abs=5, cutoff=4.5, tol=0.02):
    slopes = []
    for p in range(-max_abs, max_abs+1):
        for q in range(-max_abs, max_abs+1):
            if p == 0 or q == 0: continue
            if gcd(abs(p), abs(q)) != 1: continue
            if (p,q) in slopes or (-p,-q) in slopes: continue
            if is_blacklisted(p, q):
                continue
            slopes.append((p,q))
    
    ns = list(range(2, 61))
    occ = {n: 0 for n in ns}
    hyperbolic = []
    total = len(slopes)
    
    sys.stdout.write(f"\n{'='*60}\n")
    sys.stdout.write(f"Scanning {cusp_name} (cutoff={cutoff}, tol={tol})\n")
    sys.stdout.write(f"Blacklisted slopes: {SLOW_SLOPES}\n")
    sys.stdout.write(f"Total slopes to test: {total}\n")
    sys.stdout.write(f"{'='*60}\n")
    sys.stdout.flush()
    
    start_time = time.time()
    for idx, (p,q) in enumerate(slopes):
        slope_start = time.time()
        try:
            M = snappy.Manifold(cusp_name)
            M.dehn_fill((p,q))
            vol = float(M.volume())
            if vol < 0.5:
                elapsed = time.time() - slope_start
                sys.stdout.write(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) non-hyperbolic vol={vol:.3f} ({elapsed:.1f}s)\n")
                sys.stdout.flush()
                continue
            hits = log_n_hits(M, ns, cutoff=cutoff, tol=tol)
            hyperbolic.append((p,q,vol,len(hits)))
            for n in hits:
                occ[n] += 1
            elapsed = time.time() - slope_start
            total_elapsed = time.time() - start_time
            rate = (idx+1) / total_elapsed if total_elapsed > 0 else 0
            eta = (total - (idx+1)) / rate if rate > 0 else 0
            sys.stdout.write(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) vol={vol:.4f} hits={len(hits):2d} ({elapsed:.1f}s)  hyp={len(hyperbolic)}  eta={eta:.0f}s\n")
            sys.stdout.flush()
        except Exception as e:
            sys.stdout.write(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) ERROR: {e}\n")
            sys.stdout.flush()
            continue
    
    sys.stdout.write(f"\n{cusp_name}: {len(hyperbolic)} hyperbolic fillings out of {total}\n")
    sys.stdout.flush()
    return occ, len(hyperbolic), hyperbolic

# Run scans
sys.stdout.write("\n" + "="*60 + "\n")
sys.stdout.write("BLACKLIST SCAN: cutoff=4.5, tolerance=0.02\n")
sys.stdout.write("="*60 + "\n")
sys.stdout.flush()

m003_occ, m003_N, m003_hyp = scan_cusp('m003', max_abs=5, cutoff=4.5, tol=0.02)
m006_occ, m006_N, m006_hyp = scan_cusp('m006', max_abs=5, cutoff=4.5, tol=0.02)

# Print final tables
sys.stdout.write("\n\n" + "="*60 + "\n")
sys.stdout.write("FINAL OCCURRENCE TABLES\n")
sys.stdout.write("="*60 + "\n")
sys.stdout.flush()

sys.stdout.write(f"\nm003 (N={m003_N} hyperbolic fillings)\n")
sys.stdout.write("   n   count     %  note\n")
for n in [2,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]:
    cnt = m003_occ.get(n, 0)
    pct = 100 * cnt / m003_N if m003_N > 0 else 0
    star = "***" if cnt == m003_N else "**" if cnt > 0.8*m003_N else "*" if cnt > 0.5*m003_N else ""
    sys.stdout.write(f"  {n:3d}   {cnt:2d}/{m003_N:<2}  {pct:5.1f}%  {star}\n")
sys.stdout.flush()

sys.stdout.write(f"\nm006 (N={m006_N} hyperbolic fillings)\n")
sys.stdout.write("   n   count     %  note\n")
for n in [2,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50]:
    cnt = m006_occ.get(n, 0)
    pct = 100 * cnt / m006_N if m006_N > 0 else 0
    star = "***" if cnt == m006_N else "**" if cnt > 0.8*m006_N else "*" if cnt > 0.5*m006_N else ""
    sys.stdout.write(f"  {n:3d}   {cnt:2d}/{m006_N:<2}  {pct:5.1f}%  {star}\n")
sys.stdout.flush()

# Thresholds
thr_m003 = min([n for n in range(2,61) if m003_occ.get(n,0) > 0.8*m003_N], default=None)
thr_m006 = min([n for n in range(2,61) if m006_occ.get(n,0) > 0.8*m006_N], default=None)
sys.stdout.write(f"\nThreshold n (>80%): m003 = {thr_m003}  (log(n)/log(phi) = {np.log(thr_m003)/log_phi:.4f})\n" if thr_m003 else "\nm003: no threshold\n")
sys.stdout.write(f"Threshold n (>80%): m006 = {thr_m006}  (log(n)/log(phi) = {np.log(thr_m006)/log_phi:.4f})\n" if thr_m006 else "m006: no threshold\n")
sys.stdout.flush()

