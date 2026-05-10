import snappy
import numpy as np

M1 = snappy.OrientableClosedCensus[39]
M2 = snappy.OrientableClosedCensus[43]
M3 = snappy.OrientableClosedCensus[1]

print(f"M1: {M1.name()}, vol={float(M1.volume()):.5f}, H1={M1.homology()}")
print(f"M2: {M2.name()}, vol={float(M2.volume()):.5f}, H1={M2.homology()}")
print(f"M3: {M3.name()}, vol={float(M3.volume()):.5f}, H1={M3.homology()}")

DEPTH = 5

print(f"\nComputing ls1 (depth={DEPTH})...")
ls1 = M1.length_spectrum(DEPTH)
print(f"  done, {len(ls1)} geodesics")

print(f"Computing ls2 (depth={DEPTH})...")
ls2 = M2.length_spectrum(DEPTH)
print(f"  done, {len(ls2)} geodesics")

print(f"Computing ls3 (depth={DEPTH})...")
ls3 = M3.length_spectrum(DEPTH)
print(f"  done, {len(ls3)} geodesics")

len1 = np.array(sorted([float(abs(x.length.real)) for x in ls1]))
len2 = np.array(sorted([float(abs(x.length.real)) for x in ls2]))
len3 = np.array(sorted([float(abs(x.length.real)) for x in ls3]))

def wasserstein1(a, b):
    n = max(len(a), len(b))
    a2 = np.pad(a, (0, n - len(a)), constant_values=999.0)
    b2 = np.pad(b, (0, n - len(b)), constant_values=999.0)
    return float(np.mean(np.abs(np.sort(a2) - np.sort(b2))))

d12 = wasserstein1(len1, len2)
d13 = wasserstein1(len1, len3)
d23 = wasserstein1(len2, len3)

print(f"\nSpectral distance triangle (Wasserstein-1, depth={DEPTH}):")
print(f"  M1(census39) <-> M2/CKM(census43) : {d12:.6f}")
print(f"  M1(census39) <-> M3/PMNS(census1) : {d13:.6f}")
print(f"  M2/CKM(census43) <-> M3/PMNS(census1): {d23:.6f}")

print("\nM3 geodesic lengths (first 10):")
print([f"{l:.6f}" for l in len3[:10]])
