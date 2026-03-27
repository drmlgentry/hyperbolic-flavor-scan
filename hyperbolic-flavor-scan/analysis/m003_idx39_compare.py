import snappy
import numpy as np

# Deep comparison: m003 vs index 39
census = list(snappy.OrientableClosedCensus)
m003 = census[1]
idx39 = census[39]

print("GEODESIC SPECTRUM COMPARISON")
print("="*55)
print(f"m003  vol={float(m003.volume()):.6f}  H1={m003.homology()}")
print(f"idx39 vol={float(idx39.volume()):.6f}  H1={idx39.homology()}")
print(f"Volume ratio: {float(idx39.volume())/float(m003.volume()):.6f}")
print(f"2*vol(m003):  {2*float(m003.volume()):.6f}")
print(f"Delta:        {float(idx39.volume()) - 2*float(m003.volume()):.6f}")
print()

# Get longer spectra for both
spec_m003  = m003.length_spectrum(3.0)
spec_idx39 = idx39.length_spectrum(3.0)

m003_lens  = [(round(float(g.length.real),4),
               round(float(g.length.imag),4)) for g in spec_m003]
idx39_lens = [(round(float(g.length.real),4),
               round(float(g.length.imag),4)) for g in spec_idx39]

print(f"m003  geodesics (len<3.0): {len(m003_lens)}")
print(f"idx39 geodesics (len<3.0): {len(idx39_lens)}")

# Find shared geodesics
m003_set  = set(m003_lens)
idx39_set = set(idx39_lens)
shared    = m003_set & idx39_set
only_m003 = m003_set - idx39_set
only_39   = idx39_set - m003_set

print(f"\nShared:     {len(shared)}")
print(f"Only m003:  {len(only_m003)}")
print(f"Only idx39: {len(only_39)}")

print("\nShared geodesics:")
for g in sorted(shared):
    print(f"  len={g[0]:.4f}  phi={g[1]:.4f}")

print("\nOnly in idx39 (new geodesics):")
for g in sorted(only_39)[:10]:
    print(f"  len={g[0]:.4f}  phi={g[1]:.4f}")

# Check if idx39 volume = m003 + something meaningful
diff = float(idx39.volume()) - float(m003.volume())
print(f"\nvol(idx39) - vol(m003) = {diff:.6f}")
print(f"vol(m003)              = {float(m003.volume()):.6f}")
print(f"Ratio diff/vol(m003)   = {diff/float(m003.volume()):.6f}")

# Check Chern-Simons
print(f"\nm003  CS = {float(m003.chern_simons()):.6f}")
print(f"idx39 CS = {float(idx39.chern_simons()):.6f}")
print(f"Sum:       {float(m003.chern_simons())+float(idx39.chern_simons()):.6f}")
