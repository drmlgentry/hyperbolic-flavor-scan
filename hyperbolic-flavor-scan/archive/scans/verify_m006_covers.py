import snappy
import numpy as np

def vol(M):
    return float(M.volume())

m006 = snappy.Manifold("m006")
m003 = snappy.Manifold("m003")

print("=== POSITIVE CONTROL: idx39 -> m003 ===")
idx39 = snappy.OrientableClosedCensus[39]
print(f"idx39: vol={vol(idx39):.10f}, H1={idx39.homology()}")
print(f"m003:  vol={vol(m003):.10f},  H1={m003.homology()}")
print(f"vol ratio: {vol(idx39)/vol(m003):.10f}")

def shared_geodesics(M1, M2, max_len=3.0):
    try:
        ls1 = set(round(float(x.length.real), 6) for x in M1.length_spectrum(max_len))
        ls2 = set(round(float(x.length.real), 6) for x in M2.length_spectrum(max_len))
        return ls1 & ls2
    except Exception as e:
        return set()

shared_ctrl = shared_geodesics(idx39, m003)
print(f"Shared geodesics (len<3.0): {len(shared_ctrl)}")

print("\n=== m006 DEGREE-2 COVER CANDIDATES WITH Z/5 ===")
for name in ["m358", "s779", "m395", "s440"]:
    try:
        M = snappy.Manifold(name)
        ratio = vol(M) / vol(m006)
        shared = shared_geodesics(M, m006)
        ls_m006 = set(round(float(x.length.real),6) for x in m006.length_spectrum(3.0))
        ls_M    = set(round(float(x.length.real),6) for x in M.length_spectrum(3.0))
        subset  = ls_m006.issubset(ls_M)
        print(f"\n{name}:  vol={vol(M):.10f}  ratio={ratio:.10f}  H1={M.homology()}")
        print(f"  Shared geodesics with m006 (len<3.0): {len(shared)}")
        if shared:
            print(f"  Shared lengths: {sorted(shared)[:8]}")
        print(f"  m006 lengths subset of {name} lengths: {subset}")
    except Exception as e:
        print(f"{name}: ERROR {e}")

print("\nDone.")
