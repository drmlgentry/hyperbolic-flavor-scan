import snappy, numpy as np
from math import gcd

phi = (1+np.sqrt(5))/2
log_phi = np.log(phi)

def get_real_lengths(M, cutoff=4.5):
    spec = M.length_spectrum(cutoff)
    return [float(item.length.real()) for item in spec]

def log_n_hits(M, ns, cutoff=4.5, tol=0.02):
    try:
        lengths = get_real_lengths(M, cutoff)
        return [n for n in ns if any(abs(l-np.log(n))<tol for l in lengths)]
    except:
        return []

ns = list(range(2, 55))
results = []
seen = set()

print('Building m003 dataset...', flush=True)
for p in range(-5, 6):
    for q in range(-5, 6):
        if p==0 or q==0: continue
        if gcd(abs(p),abs(q))!=1: continue
        key = min((p,q),(-p,-q))
        if key in seen: continue
        seen.add(key)
        try:
            M = snappy.Manifold('m003')
            M.dehn_fill((p,q))
            vol = float(M.volume())
            if vol < 0.5: continue
            hits = log_n_hits(M, ns)
            results.append({'slope':(p,q),'vol':vol,'hits':hits,
                            'p_abs':abs(p),'q_abs':abs(q)})
        except:
            pass

N = len(results)
print('Loaded %d fillings' % N, flush=True)
print()

base = {n: sum(1 for r in results if n in r['hits'])/N for n in ns}

print('=== SLOPE ENCODING: p-COORDINATE ===')
print()
print('  |p|    n   N_p  N_hit  base%  p-fill%  enrich  note')
print('  ' + '-'*55)

for p_abs in [2, 3, 4, 5, 6, 7, 11, 13]:
    n = p_abs
    pf = [r for r in results if r['p_abs']==p_abs]
    if not pf:
        continue
    p_freq = sum(1 for r in pf if n in r['hits'])/len(pf)
    enr = p_freq/base[n] if base[n] > 1e-6 else 999.0
    is_luc = any(abs(phi**k + phi**(-k) - n) < 0.5 for k in range(12))
    note = '(Lucas)' if is_luc else ''
    n_hit = sum(1 for r in pf if n in r['hits'])
    print('  |p|=%2d  n=%2d  N_p=%2d  hit=%2d  %5.1f%%  %6.1f%%  %5.2fx  %s' % (
        p_abs, n, len(pf), n_hit,
        base[n]*100, p_freq*100, enr, note), flush=True)

print()
print('=== SLOPE ENCODING: q-COORDINATE ===')
print()
print('  |q|    n   N_q  N_hit  base%  q-fill%  enrich')
print('  ' + '-'*50)

for q_abs in [2, 3, 5, 7]:
    n = q_abs
    qf = [r for r in results if r['q_abs']==q_abs]
    if not qf:
        continue
    q_freq = sum(1 for r in qf if n in r['hits'])/len(qf)
    enr = q_freq/base[n] if base[n] > 1e-6 else 999.0
    n_hit = sum(1 for r in qf if n in r['hits'])
    print('  |q|=%2d  n=%2d  N_q=%2d  hit=%2d  %5.1f%%  %6.1f%%  %5.2fx' % (
        q_abs, n, len(qf), n_hit,
        base[n]*100, q_freq*100, enr), flush=True)

print()
print('=== SPARSE TIER DETAIL (n <= 17) ===')
print()
print('  n   count  base%   |p|=n enrich  note')
print('  ' + '-'*45)

for n in range(2, 18):
    cnt = sum(1 for r in results if n in r['hits'])
    freq = cnt/N
    pf = [r for r in results if r['p_abs']==n]
    if pf:
        pfreq = sum(1 for r in pf if n in r['hits'])/len(pf)
        enr = pfreq/freq if freq > 1e-6 else 999.0
        enr_str = '%4.1f%% -> %5.2fx' % (pfreq*100, enr)
    else:
        enr_str = 'no |p|=n slopes'
    k = np.log(n)/log_phi if n > 1 else 0.0
    kn = round(k*4)/4
    is_luc = abs(phi**kn + phi**(-kn) - n) < 0.5
    note = 'L_%.2f' % kn if is_luc else ''
    if freq > 0 or n in [2, 3, 7, 11]:
        print('  %2d  %2d/%2d  %5.1f%%  %-22s  %s' % (
            n, cnt, N, freq*100, enr_str, note), flush=True)
