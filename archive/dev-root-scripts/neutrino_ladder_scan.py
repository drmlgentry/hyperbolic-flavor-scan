# neutrino_ladder_scan.py

import numpy as np
import mpmath as mp
from exact_353 import build_tau_rho
from corrected_search import get_spectral_info_corrected

# Build the exact generators
tau, rho, b, c, S, G, J = build_tau_rho()

# Define the ladder word family W(m,n) = tau^m * rho * tau^n * rho^{-1}
def build_word(m, n):
    # tau^m
    if m >= 0:
        tau_m = np.linalg.matrix_power(tau, m)
    else:
        tau_m = np.linalg.matrix_power(np.linalg.inv(tau), -m)
    
    # tau^n
    if n >= 0:
        tau_n = np.linalg.matrix_power(tau, n)
    else:
        tau_n = np.linalg.matrix_power(np.linalg.inv(tau), -n)
    
    W = tau_m @ rho @ tau_n @ np.linalg.inv(rho)
    return W

# Scan range
max_pow = 15
hits = []

for m in range(-max_pow, max_pow+1):
    for n in range(-max_pow, max_pow+1):
        if m == 0 and n == 0:
            continue
        
        W = build_word(m, n)
        spec = get_spectral_info_corrected(W, high_precision=True)
        
        if spec['is_hyperbolic']:
            ell = float(spec['ell'])
            theta = spec['theta']
            # Filter for small twist and ell around 0.24
            if abs(theta) < 1e-5 and abs(ell - 0.24) < 0.1:
                hits.append((m, n, ell, theta, spec))

# Sort by closeness to ell=0.24
hits.sort(key=lambda x: abs(x[2] - 0.24))

# Print top 10
print("Top 10 candidate ladder words for neutrino sector (ell ~ 0.24, theta ~ 0):")
for m, n, ell, theta, spec in hits[:10]:
    print(f"W({m},{n}): ell={ell:.6f}, theta={theta:.3e}")

# Now, for the best candidate, do φ-quantization and minpoly analysis
if hits:
    m_best, n_best, ell_best, theta_best, spec_best = hits[0]
    print(f"\nBest candidate: W({m_best},{n_best})")
    print(f"  ell = {ell_best:.15f}")
    print(f"  theta = {theta_best:.3e}")
    
    # φ-quantization
    from phi_recognition import phi_quantization
    pq = phi_quantization(spec_best['expanding_ev'])
    print(f"  q = ell/log(phi) = {pq['q']:.15f}")
    print(f"  nearest integer = {pq['nearest_int']}")
    print(f"  residual = {pq['resid']:.3e}")
    
    # Minpoly
    from minpoly import find_minpoly
    res = find_minpoly(spec_best['expanding_ev'], prec=200, degree_max=8, maxcoeff=10**12)
    if res.found:
        print(f"  Minimal polynomial: {res.poly}")
        print(f"  Degree: {res.degree}")
        print(f"  Max error: {res.max_err:.3e}")
    else:
        print("  No minimal polynomial found.")