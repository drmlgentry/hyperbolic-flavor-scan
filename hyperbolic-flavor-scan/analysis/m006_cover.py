import snappy
import numpy as np

census = list(snappy.OrientableClosedCensus)
m006  = census[43]

print("Searching for degree-2 cover of m006...")
print(f"Target volume: 2 x {float(m006.volume()):.6f} = "
      f"{2*float(m006.volume()):.6f}")
print()

target = 2 * float(m006.volume())
m006_spec = {(round(float(g.length.real),4),
              round(float(g.length.imag),4))
             for g in m006.length_spectrum(2.0)}

candidates = []
for i, M in enumerate(census):
    try:
        vol = float(M.volume())
        if abs(vol - target) < 0.0001:
            h1 = str(M.homology())
            spec = {(round(float(g.length.real),4),
                     round(float(g.length.imag),4))
                    for g in M.length_spectrum(2.0)}
            shared = m006_spec & spec
            candidates.append((i, vol, h1, len(shared), spec))
            print(f"  [{i:5d}] vol={vol:.6f}  H1={h1:<25}  "
                  f"shared_geodesics={len(shared)}")
    except:
        pass

if not candidates:
    print("  No exact double cover found in closed census")
    print("  (may be in cusped census or requires Dehn filling)")
