import snappy
import numpy as np

M3 = snappy.OrientableClosedCensus[1]   # PMNS base
M1 = snappy.OrientableClosedCensus[39]  # degree-2 cover A
M4 = snappy.OrientableClosedCensus[40]  # degree-2 cover B (m017)

print("=== Manifold Properties ===")
for label, M in [("M3/PMNS  census[1] ", M3),
                 ("M1       census[39]", M1),
                 ("M4/m017  census[40]", M4)]:
    print(f"{label}  vol={float(M.volume()):.8f}  H1={M.homology()}")

# Compare length spectra of the two covers
print("\n=== Length Spectra Comparison (depth=4) ===")
ls1 = M1.length_spectrum(4)
ls4 = M4.length_spectrum(4)

len1 = sorted([float(abs(x.length.real)) for x in ls1])
len4 = sorted([float(abs(x.length.real)) for x in ls4])

print(f"M1 geodesic count: {len(len1)}")
print(f"M4 geodesic count: {len(len4)}")

print("\nM1 first 12 lengths:")
print([f"{l:.6f}" for l in len1[:12]])
print("M4 first 12 lengths:")
print([f"{l:.6f}" for l in len4[:12]])

tol = 1e-4
shared_14 = [(l1, l4) for l1 in len1 for l4 in len4 if abs(l1-l4) < tol]
print(f"\nShared lengths M1<->M4: {len(shared_14)}")

# Are M1 and M4 isometric? (same manifold?)
if len(shared_14) == max(len(len1), len(len4)):
    print("WARNING: M1 and M4 may be isometric!")
else:
    print("M1 and M4 are distinct covers (different length spectra)")

# Twist angles for M4 -- CP phase content
print("\n=== M4 Twist Angles (first 10) ===")
for x in ls4[:10]:
    print(f"  length={float(x.length.real):.6f}  twist={float(x.length.imag):.6f}")

# Check if M4 also has a near-zero twist geodesic like M2
twists4 = [(float(x.length.real), float(x.length.imag)) for x in ls4]
near_zero = [(l, t) for l, t in twists4 if abs(t) < 0.2]
print(f"\nNear-zero twist geodesics in M4: {near_zero}")

# Same check for M1
twists1 = [(float(x.length.real), float(x.length.imag)) for x in ls1]
near_zero1 = [(l, t) for l, t in twists1 if abs(t) < 0.2]
print(f"Near-zero twist geodesics in M1: {near_zero1}")
