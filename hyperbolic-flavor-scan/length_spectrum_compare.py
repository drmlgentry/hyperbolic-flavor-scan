import snappy

M1 = snappy.OrientableClosedCensus[39]
M2 = snappy.OrientableClosedCensus[43]

print(f"M1: {M1.name()}, vol={float(M1.volume()):.5f}, H1={M1.homology()}")
print(f"M2: {M2.name()}, vol={float(M2.volume()):.5f}, H1={M2.homology()}")

ls1 = M1.length_spectrum(4)
ls2 = M2.length_spectrum(4)

len1 = sorted([float(abs(x.length.real)) for x in ls1])
len2 = sorted([float(abs(x.length.real)) for x in ls2])

print("\nM1 geodesic lengths (first 10):")
print([f"{l:.6f}" for l in len1[:10]])

print("\nM2 geodesic lengths (first 10):")
print([f"{l:.6f}" for l in len2[:10]])

tol = 1e-4
shared = []
for l1 in len1:
    for l2 in len2:
        if abs(l1 - l2) < tol:
            shared.append((l1, l2, abs(l1 - l2)))

print("\nShared geodesic lengths (within tol=1e-4):")
if shared:
    for a, b, diff in shared:
        print(f"  M1={a:.6f}  M2={b:.6f}  diff={diff:.2e}")
else:
    print("  None found at this depth -- try increasing length_spectrum depth")

print("\nM1 (length, twist_angle) pairs:")
for x in ls1[:10]:
    print(f"  length={float(x.length.real):.6f}  twist={float(x.length.imag):.6f}")

print("\nM2 (length, twist_angle) pairs:")
for x in ls2[:10]:
    print(f"  length={float(x.length.real):.6f}  twist={float(x.length.imag):.6f}")
