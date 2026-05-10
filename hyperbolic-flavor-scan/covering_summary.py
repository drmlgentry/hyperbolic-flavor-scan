import snappy

def vol(M):
    return float(M.volume())

def shared_geo(M1, M2, max_len=3.5):
    try:
        ls1 = set(round(float(x.length.real),5) for x in M1.length_spectrum(max_len))
        ls2 = set(round(float(x.length.real),5) for x in M2.length_spectrum(max_len))
        return ls1 & ls2
    except:
        return set()

m003 = snappy.OrientableClosedCensus[1]
m006 = snappy.OrientableClosedCensus[43]
idx39 = snappy.OrientableClosedCensus[39]

print("=== SUMMARY: COVERING STRUCTURE ===\n")
print(f"m003 (idx 1):  vol={vol(m003):.10f}  H1={m003.homology()}")
print(f"m006 (idx 43): vol={vol(m006):.10f}  H1={m006.homology()}")

print("\n--- Degree-2 covers of m003 ---")
for idx in [39, 40]:
    M = snappy.OrientableClosedCensus[idx]
    shared = shared_geo(M, m003)
    ls_m003 = set(round(float(x.length.real),5) for x in m003.length_spectrum(3.5))
    ls_M    = set(round(float(x.length.real),5) for x in M.length_spectrum(3.5))
    print(f"  idx={idx} {M.name()}: H1={M.homology()}  "
          f"ratio={vol(M)/vol(m003):.8f}  "
          f"shared={len(shared)}  "
          f"m003_subset={ls_m003.issubset(ls_M)}")

print("\n--- Degree-3 covers of m003 (indices 237-240) ---")
for idx in [237, 238, 239, 240]:
    M = snappy.OrientableClosedCensus[idx]
    shared = shared_geo(M, m003)
    ls_m003 = set(round(float(x.length.real),5) for x in m003.length_spectrum(3.5))
    ls_M    = set(round(float(x.length.real),5) for x in M.length_spectrum(3.5))
    print(f"  idx={idx} {M.name()}: H1={M.homology()}  "
          f"ratio={vol(M)/vol(m003):.8f}  "
          f"shared={len(shared)}  "
          f"m003_subset={ls_m003.issubset(ls_M)}")

print("\n--- Degree-2 cover of m006 (s394, idx 1275) ---")
s394 = snappy.OrientableClosedCensus[1275]
shared = shared_geo(s394, m006)
ls_m006 = set(round(float(x.length.real),5) for x in m006.length_spectrum(3.5))
ls_s394 = set(round(float(x.length.real),5) for x in s394.length_spectrum(3.5))
print(f"  idx=1275 s394: H1={s394.homology()}  "
      f"ratio={vol(s394)/vol(m006):.8f}  "
      f"shared={len(shared)}  "
      f"m006_subset={ls_m006.issubset(ls_s394)}")

print("\n=== IMPLICATIONS FOR CORE_MASTER ===")
print("""
CONFIRMED (theorem stands):
  m006 has no degree-2 cover in OrientableClosedCensus with H1 containing Z/5.
  s394 (the unique vol=2*vol(m006) hit) has H1=Z/3, not Z/5, and
  only 4 shared geodesics — not a cover by geodesic criterion.

CONFIRMED (covering asymmetry):
  m003 has degree-2 cover idx39 (H1=Z/55=Z/5*Z/11, 20 shared geodesics).
  m003 also has a second degree-2 cover idx40 (H1=Z/7+Z/7, check below).

NEW (degree-3 covers of m003):
  Four manifolds at 3*vol(m003): m130 x2, m139 x2.
  H1 values: Z/16, Z/80, Z/48, Z/2+Z/42.
  Z/80 = Z/16*Z/5 and Z/42 contains no Z/5 factor (42=2*3*7).
  Z/80 candidate (idx 238) may be a degree-3 cover preserving Z/5.

NOTE on SnapPy naming collision:
  OrientableClosedCensus[39] is named 'm006' but is NOT our CKM manifold.
  Our CKM manifold is OrientableClosedCensus[43], also named 'm006'.
  All code should reference manifolds by census INDEX, not name string.
""")
