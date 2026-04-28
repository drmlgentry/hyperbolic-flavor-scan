import snappy
import numpy as np
import time
from math import gcd

phi = (1 + np.sqrt(5)) / 2
log_phi = np.log(phi)

# Blacklist slopes that are known to be slow or problematic
SLOW_SLOPES = [
    (-5, 4), (5, -4),   # from earlier run – hung for a long time
    (-4, 5), (4, -5),   # also slow
    # Add more as you encounter them
]

def is_blacklisted(p, q):
    return (p, q) in SLOW_SLOPES or (-p, -q) in SLOW_SLOPES

def get_real_lengths(M, cutoff=4.5, timeout=60):
    """Get lengths with a simple timeout (approximate)."""
    start = time.time()
    try:
        spec = M.length_spectrum(cutoff)
        elapsed = time.time() - start
        if elapsed > timeout:
            print(f"      WARNING: length_spectrum took {elapsed:.1f}s (>{timeout}s)")
        return [float(item.length.real()) for item in spec]
    except Exception as e:
        print(f"      ERROR in length_spectrum: {e}")
        return []

def log_n_hits(M, ns, cutoff=4.5, tol=0.02):
    lengths = get_real_lengths(M, cutoff)
    if not lengths:
        return []
    return [n for n in ns if any(abs(l - np.log(n)) < tol for l in lengths)]

def scan_fillings(max_abs=5, cutoff=4.5, tol=0.02):
    slopes = []
    seen = set()
    for p in range(-max_abs, max_abs + 1):
        for q in range(-max_abs, max_abs + 1):
            if p == 0 or q == 0:
                continue
            if gcd(abs(p), abs(q)) != 1:
                continue
            key = min((p, q), (-p, -q))
            if key in seen:
                continue
            seen.add(key)
            slopes.append((p, q))
    
    ns = list(range(2, 55))
    results = []
    hyperbolic = []
    total = len(slopes)
    
    print(f"Total slopes to test: {total}")
    print(f"Blacklisted slopes: {SLOW_SLOPES}")
    print()
    
    for idx, (p, q) in enumerate(slopes):
        if is_blacklisted(p, q):
            print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) SKIPPED (blacklisted)")
            continue
        
        slope_start = time.time()
        try:
            M = snappy.Manifold('m003')
            M.dehn_fill((p, q))
            vol = float(M.volume())
            if vol < 0.5:
                elapsed = time.time() - slope_start
                print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) non-hyperbolic vol={vol:.3f} ({elapsed:.1f}s)")
                continue
            
            hits = log_n_hits(M, ns, cutoff=cutoff, tol=tol)
            hyperbolic.append((p, q, vol, len(hits)))
            for n in hits:
                results.append((p, q, n))
            elapsed = time.time() - slope_start
            print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) vol={vol:.4f} hits={len(hits):2d} ({elapsed:.1f}s)  hyp={len(hyperbolic)}")
        except Exception as e:
            print(f"[{idx+1:3d}/{total}] ({p:2d},{q:2d}) ERROR: {e}")
            continue
    
    return results, hyperbolic

# Run scan
print("=" * 60)
print("SLOPE ENCODING SCAN: m003")
print("cutoff=4.5, tolerance=0.02")
print("=" * 60)

results, hyperbolic = scan_fillings(max_abs=5, cutoff=4.5, tol=0.02)

# Analyse slope encoding
print("\n" + "=" * 60)
print("SLOPE ENCODING ANALYSIS")
print("=" * 60)

# Baseline occurrence for each n
from collections import defaultdict
total_hits = defaultdict(int)
for p, q, n in results:
    total_hits[n] += 1
N = len(hyperbolic)

print(f"\nTotal hyperbolic fillings: {N}")
print("\nBaseline occurrence (all fillings):")
for n in [2, 3, 5, 7, 11, 13, 17, 19]:
    cnt = total_hits.get(n, 0)
    pct = 100 * cnt / N if N > 0 else 0
    print(f"  n={n:2d}: {cnt:2d}/{N} ({pct:5.1f}%)")

# Enrichment for |p| = n
print("\nEnrichment for fillings with |p| = n:")
for target_n in [2, 3, 5, 7, 11, 13]:
    # Count fillings with |p| = target_n
    p_fillings = [(p, q, v, h) for (p, q, v, h) in hyperbolic if abs(p) == target_n]
    Np = len(p_fillings)
    if Np == 0:
        print(f"  |p|={target_n}: no fillings")
        continue
    
    # Count hits for that target n in those fillings
    hits_p = 0
    for p, q, n in results:
        if abs(p) == target_n and n == target_n:
            hits_p += 1
    
    # Baseline hits for the same n across all fillings
    hits_all = total_hits.get(target_n, 0)
    
    baseline_pct = 100 * hits_all / N if N > 0 else 0
    p_pct = 100 * hits_p / Np if Np > 0 else 0
    enrichment = p_pct / baseline_pct if baseline_pct > 0 else 0
    
    print(f"  |p|={target_n}: Np={Np:2d}, hits_p={hits_p:2d}, p_pct={p_pct:5.1f}%, baseline={baseline_pct:5.1f}%, enrichment={enrichment:.2f}x")

# Check for |p| = 13 (predicted no enrichment)
print("\nPrediction: |p| = 13 should show no enrichment (baseline ~0)")

print("\n" + "=" * 60)
print("SCAN COMPLETE")
print("=" * 60)

