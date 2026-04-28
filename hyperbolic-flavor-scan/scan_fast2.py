import snappy
import numpy as np
from math import gcd
import time

phi = (1+np.sqrt(5))/2
log_phi = np.log(phi)

def get_real_lengths(M, cutoff=5.0):
    spec = M.length_spectrum(cutoff)
    return [float(item.length.real()) for item in spec]

def log_n_hits(M, ns, cutoff=5.0, tol=0.02):
    try:
        lengths = get_real_lengths(M, cutoff)
        return [n for n in ns
                if any(abs(l - np.log(n)) < tol for l in lengths)]
    except:
        return []

def scan_cusp(cusp_name, max_abs=5):
    slopes = []
    for p in range(-max_abs, max_abs+1):
        for q in range(-max_abs, max_abs+1):
            if p == 0 or q == 0: continue
            if gcd(abs(p), abs(q)) != 1: continue
            if (p,q) in slopes or (-p,-q) in slopes: continue
            slopes.append((p,q))

    ns = list(range(2, 61))
    results = {}
    total = len(slopes)

    print(f"\n=== Scanning {cusp_name} ({total} slopes, cutoff=5.0) ===",
          flush=True)
    start_time = time.time()

    for idx, (p,q) in enumerate(slopes):
        try:
            M = snappy.Manifold(cusp_name)
            M.dehn_fill((p,q))
            vol = float(M.volume())
            if vol < 0.5:
                print(f"  [{idx+1:>3}/{total}] ({p:>3},{q:>3}): "
                      f"non-hyperbolic (vol={vol:.3f})", flush=True)
                continue
            hits = log_n_hits(M, ns)
            results[(p,q)] = {'hits': hits, 'vol': vol,
                              'p_abs': abs(p), 'q_abs': abs(q)}
            print(f"  [{idx+1:>3}/{total}] ({p:>3},{q:>3}): "
                  f"vol={vol:.4f}  {len(hits):>3} hits", flush=True)
        except Exception as e:
            print(f"  [{idx+1:>3}/{total}] ({p:>3},{q:>3}): error {e}",
                  flush=True)

        if (idx+1) % 10 == 0:
            elapsed = time.time() - start_time
            rate = (idx+1)/elapsed
            eta = (total-(idx+1))/rate
            print(f"    --> progress: {idx+1}/{total} "
                  f"({100*(idx+1)/total:.0f}%)  "
                  f"elapsed={elapsed:.0f}s  ETA={eta:.0f}s",
                  flush=True)

    return results, ns

# Run
print("="*60, flush=True)
print("FAST SCAN v2: cutoff=5.0, |p|,|q|<=5", flush=True)
print("="*60, flush=True)

m003_res, ns = scan_cusp('m003')
m006_res, _  = scan_cusp('m006')

N3 = len(m003_res)
N6 = len(m006_res)

print(f"\n{'='*60}", flush=True)
print(f"DONE: m003={N3} fillings  m006={N6} fillings", flush=True)
print(f"{'='*60}", flush=True)

# Occurrence frequencies
for cusp, res, N in [('m003', m003_res, N3), ('m006', m006_res, N6)]:
    print(f"\n=== {cusp} OCCURRENCE FREQUENCY ({N} fillings) ===",
          flush=True)
    print(f"  {'n':>4} {'count':>7} {'%':>6}  note", flush=True)
    print("  " + "-"*38, flush=True)
    for n in ns:
        cnt = sum(1 for r in res.values() if n in r['hits'])
        pct = 100*cnt/N
        k = np.log(n)/log_phi
        k4 = round(k*4)/4
        is_luc = abs(phi**k4 + phi**(-k4) - n) < 0.5
        if pct > 5 or is_luc:
            note = f"L_{k4:.2f}" if is_luc else ""
            star = "***" if cnt==N else "**" if pct>80 else "*" if pct>50 else ""
            print(f"  {n:>4} {cnt:>4}/{N:<3} {pct:>5.1f}%  "
                  f"{star:<4} {note}", flush=True)

# Slope encoding
print(f"\n=== SLOPE ENCODING m003 ({N3} fillings) ===", flush=True)
base = {n: sum(1 for r in m003_res.values() if n in r['hits'])/N3
        for n in ns}
print(f"  {'|p|':>5} {'n':>4} {'N_p':>5} {'base%':>7} "
      f"{'p%':>7} {'enrich':>8}  note", flush=True)
print("  "+"-"*52, flush=True)
for p_abs in [2,3,5,6,7,11,13]:
    n = p_abs
    pf = [r for r in m003_res.values() if r['p_abs']==p_abs]
    if len(pf) < 2: continue
    pfreq = sum(1 for r in pf if n in r['hits'])/len(pf)
    enr = pfreq/base[n] if base[n]>1e-6 else 999
    is_luc = any(abs(phi**k+phi**(-k)-n)<0.5 for k in range(12))
    note = "(Lucas)" if is_luc else ""
    print(f"  {p_abs:>5} {n:>4} {len(pf):>5} {base[n]*100:>6.1f}% "
          f"{pfreq*100:>6.1f}% {enr:>7.2f}x  {note}", flush=True)

# Phase transition
print(f"\n=== PHASE TRANSITION ===", flush=True)
for cusp, res, N in [('m003', m003_res, N3), ('m006', m006_res, N6)]:
    freq = {n: sum(1 for r in res.values() if n in r['hits'])/N
            for n in ns}
    univ  = [n for n in ns if freq[n] >= 0.80]
    mid   = [n for n in ns if 0.40 <= freq[n] < 0.80]
    sparse= [n for n in ns if 0 < freq[n] < 0.40]
    lk    = [i/4 for i in range(40)]
    lu = [n for n in univ   if any(abs(phi**k+phi**(-k)-n)<0.5 for k in lk)]
    ls = [n for n in sparse if any(abs(phi**k+phi**(-k)-n)<0.5 for k in lk)]
    print(f"\n{cusp} ({N} fillings):", flush=True)
    print(f"  Universal (>=80%): first={min(univ) if univ else 'none'}"
          f"  count={len(univ)}", flush=True)
    print(f"  Mid (40-80%):  {mid[:10]}", flush=True)
    print(f"  Sparse (<40%): {sparse[:10]}", flush=True)
    print(f"  Lucas universal: {lu[:10]}", flush=True)
    print(f"  Lucas sparse:    {ls[:8]}", flush=True)
