import snappy
import numpy as np

PHI = (1+5**0.5)/2
vol_m003 = 0.9814
vol_m006 = 2.0289

targets = [
    (1,   "m003"),
    (43,  "m006"),
    (39,  "2xm003 candidate"),
    (431, "m006*phi candidate"),
    (246, "m003+m006 candidate (close)"),
    (254, "m003+m006 candidate (close)"),
]

census = list(snappy.OrientableClosedCensus)

for idx, label in targets:
    M = census[idx]
    print(f"\n[{idx}] {label}")
    print(f"  Volume:  {float(M.volume()):.6f}")
    print(f"  H1:      {M.homology()}")
    try:
        print(f"  Chern-Simons: {M.chern_simons():.6f}")
    except: pass
    try:
        spec = M.length_spectrum(1.5)
        print(f"  Short geodesics (len<1.5): {len(spec)}")
        for g in spec[:5]:
            print(f"    len={float(g.length.real):.4f}  "
                  f"phi={float(g.length.imag):.4f}")
    except Exception as e:
        print(f"  Spectrum error: {e}")
