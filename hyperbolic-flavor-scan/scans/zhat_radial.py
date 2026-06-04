#!/usr/bin/env python3
import numpy as np
import mpmath as mp
mp.mp.dps = 30
def eta_product(q, N_max=200):
    q = mp.mpc(q)
    f = q
    for n in range(1, N_max+1):
        f *= (1 - q**n)**2 * (1 - q**(11*n))**2
    return f
def radial_limit(r, N_max=200):
    q = mp.exp(mp.mpc(0, 2*np.pi/r))
    return eta_product(q, N_max)
def test_wrt(r):
    val = radial_limit(r)
    abs2 = float(abs(val)**2)
    target = 13.0 / r
    return abs2, target, (abs2/target - 1.0)*100
print("r\t|Z_hat(r)|^2\t13/r\t\tdiff %")
for r in range(1,46,2):
    if r%5==0: continue
    try:
        a,t,e = test_wrt(r)
        if a>1e-10: print(f"{r}\t{a:.6f}\t{t:.6f}\t{e:.4f}%")
    except: pass
# NOTE: Naive eta product collapses at roots of unity (exact zeros in product).
# The Z-hat conjecture requires a different computational approach (GPPV formula).
