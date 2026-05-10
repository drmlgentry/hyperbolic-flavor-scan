import snappy
import numpy as np

M1 = snappy.OrientableClosedCensus[39]
M3 = snappy.OrientableClosedCensus[1]

v1 = float(M1.volume())
v3 = float(M3.volume())

print(f"M1 volume: {v1:.8f}")
print(f"M3 volume: {v3:.8f}")
print(f"Ratio M1/M3: {v1/v3:.8f}")
print(f"Ratio M3/M1: {v3/v1:.8f}")

# Compare length spectra directly for shared lengths
ls1 = M1.length_spectrum(5)
ls3 = M3.length_spectrum(5)

len1 = sorted([float(abs(x.length.real)) for x in ls1])
len3 = sorted([float(abs(x.length.real)) for x in ls3])

tol = 1e-4
shared = []
for l1 in len1:
    for l3 in len3:
        if abs(l1 - l3) < tol:
            shared.append((l1, l3, abs(l1-l3)))

print(f"\nShared geodesic lengths M1 <-> M3 (tol=1e-4):")
print(f"  Count: {len(shared)} of {len(len1)} M1 geodesics")
for a, b, d in shared[:20]:
    print(f"  M1={a:.6f}  M3={b:.6f}  diff={d:.2e}")

# Check if M3 lengths appear scaled by 2 in M1
print(f"\nChecking if M3 lengths * 2 appear in M1:")
for l3 in len3[:15]:
    target = 2.0 * l3
    matches = [l1 for l1 in len1 if abs(l1 - target) < tol]
    if matches:
        print(f"  2 * M3({l3:.6f}) = {target:.6f} found in M1: {matches}")
