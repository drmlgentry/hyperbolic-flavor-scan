import snappy
import numpy as np

# For arithmetic manifolds, the CS invariant satisfies:
# CS(M) = eta(M)/2 mod 1
# where eta is the eta invariant of the Dirac operator
# For m003: CS = 1/4, so eta(m003) = 1/2 mod 2

# The Dedekind sum relation (Atiyah-Patodi-Singer):
# CS(M_p/q) = CS(M) + s(p,q) mod 1
# where s(p,q) is the Dedekind sum and M_p/q is Dehn filling

# Compute the Dehn filling that gives OrientableClosedCensus[1]
M_cusp = snappy.Manifold("m003")
print(f"m003 cusped: {M_cusp}")
print(f"CS(cusped) = {float(M_cusp.chern_simons()):.8f}")
print(f"Dehn fillings: {M_cusp.dehn_filling_matrix()}" 
      if hasattr(M_cusp,"dehn_filling_matrix") else "")

# Try to identify the Dehn surgery coefficients
M_closed = snappy.OrientableClosedCensus[1]
print(f"\nm003 closed vol = {float(M_closed.volume()):.8f}")

# The surgery description
try:
    surgery = M_cusp.hyperbolic_dehn_surgery((1,0))
    print(f"(1,0) filling vol = {float(surgery.volume()):.8f}")
except Exception as e:
    print(f"Surgery test: {e}")

# Scan Dehn fillings to find the closed manifold
print("\nScanning Dehn fillings of m003 cusped to find OrientableClosedCensus[1]:")
target_vol = 0.98136883
for p in range(-10, 11):
    for q in range(-10, 11):
        if q == 0: continue
        from math import gcd
        if abs(gcd(p,q)) != 1: continue
        try:
            M_fill = M_cusp.dehn_fill((p,q))
            v = float(M_fill.volume())
            if abs(v - target_vol) < 0.001:
                cs_fill = float(M_fill.chern_simons()) if M_fill.is_orientable() else None
                print(f"  ({p:3d},{q:3d}) vol={v:.6f} CS={cs_fill:.6f if cs_fill else 'N/A'}")
        except: pass
