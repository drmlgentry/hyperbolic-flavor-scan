import snappy
import time
from math import gcd

# Test a few specific slopes that were fast in your earlier run
test_slopes = [
    (-5, -3),
    (-5, -2), 
    (-5, -1),
    (-5, 1),
    (-5, 2),
    (-5, 3),
    (-4, -3),
    (-4, -1),
    (-3, -5),
    (-3, -4),
    (-3, -2),
    (-2, -5),
    (-2, -3),
    (-2, -1),
    (-1, -5),
    (-1, -4),
    (-1, -3),
]

print("Testing individual slopes for compute time...")
print("Slope\t\tTime(s)\tVol\tHits")
print("-" * 50)

for p, q in test_slopes:
    start = time.time()
    try:
        M = snappy.Manifold('m003')
        M.dehn_fill((p, q))
        vol = float(M.volume())
        if vol < 0.5:
            elapsed = time.time() - start
            print(f"({p:2d},{q:2d})\t\t{elapsed:.1f}s\t{vol:.3f}\tnon-hyp")
            continue
        spec = M.length_spectrum(4.5)
        lengths = len(spec)
        elapsed = time.time() - start
        print(f"({p:2d},{q:2d})\t\t{elapsed:.1f}s\t{vol:.3f}\t{lengths}")
    except Exception as e:
        elapsed = time.time() - start
        print(f"({p:2d},{q:2d})\t\t{elapsed:.1f}s\tERROR: {e}")

