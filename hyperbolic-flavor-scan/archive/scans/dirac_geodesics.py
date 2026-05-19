import snappy
import numpy as np
import math

PHI = (1 + math.sqrt(5)) / 2
LOG_PHI = math.log(PHI)

for idx, name in [(1, "m003"), (43, "m006")]:
    M = snappy.OrientableClosedCensus[idx]
    print()
    print("="*50)
    print(f"{name} = OrientableClosedCensus[{idx}]")
    print(f"Volume: {M.volume():.10f}")
    spec = M.length_spectrum(cutoff=3.0)
    print(f"Geodesics up to length 3.0: {len(spec)}")
    print()
    print(f"  {chr(39)}length{chr(39):>12}  {chr(39)}mult{chr(39):>6}  {chr(39)}twist{chr(39):>10}  {chr(39)}l/log(phi){chr(39):>12}")
    print("-"*55)
    for g in spec[:30]:
        ell = float(g.length.real)
        twist = float(g.length.imag)
        mult = g.multiplicity
        ratio = ell / LOG_PHI
        q = round(ratio * 4) / 4
        res = abs(ratio - q)
        flag = "***" if res < 0.02 else "**" if res < 0.05 else "*" if res < 0.10 else ""
        print(f"  {ell:12.6f}  {mult:6d}  {twist:10.6f}  {ratio:12.4f}  q={q:+.2f}  res={res:.4f} {flag}")
