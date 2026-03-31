import snappy
import numpy as np

def vol(M):
    return float(M.volume())

def shared_geo(M1, M2, max_len=3.5):
    try:
        ls1 = set(round(float(x.length.real), 5) for x in M1.length_spectrum(max_len))
        ls2 = set(round(float(x.length.real), 5) for x in M2.length_spectrum(max_len))
        return ls1 & ls2
    except:
        return set()

# Correct anchors — closed census only
m003 = snappy.OrientableClosedCensus[1]   # vol=0.9814, H1=Z/5, PMNS
m006 = snappy.OrientableClosedCensus[43]  # vol=2.0289, H1=Z/5, CKM
idx39 = snappy.OrientableClosedCensus[39] # vol=1.9627, H1=Z/55, known m003 cover

VOL_m003 = vol(m003)
VOL_m006 = vol(m006)
PHI = (1 + 5**0.5) / 2

print("=== ANCHORS ===")
print(f"m003  (idx 1):  vol={VOL_m003:.10f}  H1={m003.homology()}")
print(f"m006  (idx 43): vol={VOL_m006:.10f}  H1={m006.homology()}")
print(f"idx39 (idx 39): vol={vol(idx39):.10f}  H1={idx39.homology()}")
print(f"idx39/m003 ratio: {vol(idx39)/VOL_m003:.10f}")

print("\n=== POSITIVE CONTROL: idx39 -> m003 ===")
shared = shared_geo(idx39, m003)
ls_m003 = set(round(float(x.length.real),5) for x in m003.length_spectrum(3.5))
ls_idx39 = set(round(float(x.length.real),5) for x in idx39.length_spectrum(3.5))
print(f"Shared geodesics (len<3.5): {len(shared)}")
print(f"m003 lengths subset of idx39: {ls_m003.issubset(ls_idx39)}")
if shared:
    print(f"Shared lengths: {sorted(shared)[:10]}")

print("\n=== TARGET VOLUMES (CORRECTED) ===")
targets = {
    "2*vol(m003)  [idx39, known]":  2 * VOL_m003,
    "2*vol(m006)":                  2 * VOL_m006,
    "3*vol(m003)":                  3 * VOL_m003,
    "vol(m003)+vol(m006)":          VOL_m003 + VOL_m006,
    "phi*vol(m003)":                PHI * VOL_m003,
    "phi*vol(m006)":                PHI * VOL_m006,
}
for k,v2 in targets.items():
    print(f"  {k:40s} = {v2:.8f}")

print("\n=== CENSUS SCAN WITH CORRECTED VOLUMES ===")
census = snappy.OrientableClosedCensus
tol = 0.0001
hits = []
for i, M in enumerate(census):
    v2 = vol(M)
    for label, target in targets.items():
        if abs(v2 - target) < tol:
            h1 = str(M.homology())
            hits.append((label, i, M.name(), v2, h1))

# Print sorted by target
for label in targets:
    group = [(i,n,v2,h1) for l,i,n,v2,h1 in hits if l==label]
    if group:
        print(f"\n  {label}:")
        for i,n,v2,h1 in group:
            print(f"    idx={i:5d}  name={n:8s}  vol={v2:.10f}  H1={h1}")
    else:
        print(f"\n  {label}: NO HITS")

print("\n=== VERIFY m006 COVER CANDIDATES ===")
m006_covers = [(i,n,v2,h1) for l,i,n,v2,h1 in hits if "2*vol(m006)" in l]
for i,n,v2,h1 in m006_covers:
    M = census[i]
    shared = shared_geo(M, m006)
    ls_m006 = set(round(float(x.length.real),5) for x in m006.length_spectrum(3.5))
    ls_M    = set(round(float(x.length.real),5) for x in M.length_spectrum(3.5))
    subset  = ls_m006.issubset(ls_M)
    z5      = "5" in h1
    print(f"\n  idx={i} {n}: H1={h1}  has_Z5={z5}")
    print(f"    shared geodesics={len(shared)}  m006_lengths_subset={subset}")
    if shared:
        print(f"    shared: {sorted(shared)[:6]}")

print("\nDone.")
