import snappy, numpy as np

PHI = (1+5**0.5)/2

print("=== Full spectral estimate: first 20 geodesics ===")
print("Using lambda_approx = 1/4 + (2*pi*n/l)^2 for n=1")
print("This is an approximation, not the true spectrum.")
print()

# Null test: how many hits expected by chance?
# Each estimate uniform in log_phi, so P(resid < 0.05) = 0.05/0.125 = 0.40
# With N estimates, expected hits = 0.40 * N
# Need significantly more than 40% to claim structure

for idx, name in [(1,'PMNS'),(43,'CKM')]:
    M = snappy.OrientableClosedCensus[idx]
    ls = M.length_spectrum(4.0)
    geos = list(ls)
    
    hits = 0
    total = 0
    print(f"{name}:")
    for g in geos[:20]:
        l = float(g['length'].real())
        n = 1
        lam = 0.25 + (2*np.pi*n/l)**2
        log_lam = np.log(lam)/np.log(PHI)
        nearest_q = round(log_lam/0.25)*0.25
        residual = abs(log_lam - nearest_q)
        total += 1
        if residual < 0.05:
            hits += 1
            print(f"  HIT: l={l:.4f} lambda={lam:.2f} resid={residual:.4f}")
    
    hit_rate = hits/total
    expected = 0.40
    print(f"  Hit rate: {hits}/{total} = {hit_rate:.2f} "
          f"(expected by chance: {expected:.2f})")
    print(f"  Excess above chance: {hit_rate - expected:+.2f}")
    print()
