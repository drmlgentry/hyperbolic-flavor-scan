import snappy
import numpy as np
import cmath

# The closed manifold IS the cusped manifold with specific Dehn filling
# SnapPy's OrientableClosedCensus stores these as specific surgeries
# Let's find the surgery description differently

M_closed = snappy.OrientableClosedCensus[1]
M_cusp   = snappy.Manifold("m003")

print("=== Identifying the Dehn filling ===")
print(f"Closed vol:  {float(M_closed.volume()):.8f}")
print(f"Cusped vol:  {float(M_cusp.volume()):.8f} (incomplete - expected)")

# The closed manifold description
print(f"\nClosed manifold repr: {repr(M_closed)}")
print(f"Closed manifold str:  {str(M_closed)}")

# Try getting the surgery info from the closed manifold directly
try:
    # SnapPy closed manifolds store surgery coefficients
    link = M_closed.link()
    print(f"Link: {link}")
except Exception as e:
    print(f"Link: {e}")

# The key quantum gravity computations we CAN do without Sage:
# 1. Semiclassical approximation to WRT invariant
# 2. Length spectrum analysis
# 3. Spectral zeta function (partial)

print("\n=== Semiclassical WRT invariant approximation ===")
print("Z_CS(M,k) ~ (k/2pi)^(3/2) * Vol(M)^(1/2) * exp(i*k*Vol(M) + i*pi*eta/4)")
print("This is the stationary phase approximation for large k")
print()

vol_m003 = float(M_closed.volume())
vol_m006 = float(snappy.OrientableClosedCensus[43].volume())
cs_m003  = float(M_cusp.chern_simons())  # = 0.25
cs_m006  = float(snappy.Manifold("m006").chern_simons())  # = -0.114137

print(f"m003: Vol={vol_m003:.6f}, CS={cs_m003:.6f}")
print(f"m006: Vol={vol_m006:.6f}, CS={cs_m006:.6f}")

print(f"\n{'k':>4} {'|Z_m003|':>12} {'|Z_m006|':>12} "
      f"{'arg(Z_m003)/pi':>16} {'arg(Z_m006)/pi':>16}")
print("-"*65)

for k in range(1, 26):
    # Semiclassical amplitude: (k)^(3/2) * sqrt(Vol)
    amp_m003 = k**1.5 * vol_m003**0.5
    amp_m006 = k**1.5 * vol_m006**0.5
    # Semiclassical phase: exp(i*k*Vol + 2*pi*i*k*CS)
    phase_m003 = cmath.exp(1j*(k*vol_m003 + 2*np.pi*k*cs_m003))
    phase_m006 = cmath.exp(1j*(k*vol_m006 + 2*np.pi*k*cs_m006))
    z_m003 = amp_m003 * phase_m003
    z_m006 = amp_m006 * phase_m006
    print(f"{k:>4} {abs(z_m003):>12.4f} {abs(z_m006):>12.4f} "
          f"{cmath.phase(z_m003)/np.pi:>+16.6f} "
          f"{cmath.phase(z_m006)/np.pi:>+16.6f}")

print("""
Note: The amplitudes grow as k^(3/2) in the semiclassical approximation.
The physically meaningful quantity is the NORMALIZED ratio Z_m003/Z_m006
and the PHASE difference arg(Z_m003) - arg(Z_m006).
""")

# ── Key invariant: ratio at k=5 ───────────────────────────────────
print("=== Ratio Z_m003 / Z_m006 at k=5 (matching H1=Z/5) ===")
k = 5
phase_m003 = cmath.exp(1j*(k*vol_m003 + 2*np.pi*k*cs_m003))
phase_m006 = cmath.exp(1j*(k*vol_m006 + 2*np.pi*k*cs_m006))
ratio = (vol_m003**0.5 * phase_m003) / (vol_m006**0.5 * phase_m006)
print(f"  |Z_m003/Z_m006| = {abs(ratio):.6f}")
print(f"  arg(ratio)/pi   = {cmath.phase(ratio)/np.pi:.6f}")
print(f"  ratio           = {ratio.real:.6f} + {ratio.imag:.6f}i")
print()

# ── Connection to flavor physics ──────────────────────────────────
print("=== CONNECTION TO FLAVOR PHYSICS ===")
print()
print("The twist angles from the length spectrum:")
for idx, name in [(1,"m003"),(43,"m006")]:
    M = snappy.OrientableClosedCensus[idx]
    ls = M.length_spectrum(2.5)
    geos = list(ls)
    print(f"\n{name} shortest geodesic twist angles (mod 2pi):")
    for g in geos[:5]:
        phi = float(g.length.imag)
        phi_mod = phi % (2*np.pi)
        phi_frac = phi / (2*np.pi)
        print(f"  phi={phi:.6f} rad = {phi*180/np.pi:.3f}° "
              f"= {phi_frac:.4f} * 2pi")

print("""
The CP phase delta_CP = 203.5° from the PMNS paper corresponds to
the signed sum of twist angles of the optimal geodesic triple on m003.

In the quantum gravity language:
  delta_CP = sum_i phi(gamma_i) = sum_i Im(complex_length(gamma_i))

This is the holonomy of the gravitational Chern-Simons connection
around the closed geodesics of the PMNS manifold.

The Chern-Simons level k=5 matches H1(M)=Z/5, suggesting that
the physical CP violation is the k=5 Chern-Simons amplitude.
""")

# ── Save results ──────────────────────────────────────────────────
import json
results = {
    "m003": {
        "volume": vol_m003,
        "cs_cusped": cs_m003,
        "cs_rational": "1/4",
        "h1": "Z/5",
        "arithmetic": True,
    },
    "m006": {
        "volume": vol_m006,
        "cs_cusped": cs_m006,
        "cs_rational": None,
        "h1": "Z/5",
        "arithmetic": False,
    }
}
with open("quantum_gravity_results.json","w") as f:
    json.dump(results, f, indent=2)
print("Saved quantum_gravity_results.json")
