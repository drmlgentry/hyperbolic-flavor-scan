import snappy
import numpy as np

PHI = (1 + 5**0.5) / 2
TARGET_VOL = 0.9814 + 2.0289  # 3.0103

print("Searching OrientableClosedCensus for parent A5 manifold...")
print(f"Target volume: {TARGET_VOL:.4f}")
print(f"(vol(m003) + vol(m006) = 0.9814 + 2.0289)")
print()

candidates = []
census = snappy.OrientableClosedCensus

for i, M in enumerate(census):
    try:
        vol = float(M.volume())
        h1  = str(M.homology())
        
        # Check Z/5 torsion
        has_z5 = "Z/5" in h1 or "Z5" in h1
        
        # Check volume proximity to target
        # Also check multiples and simple combinations
        targets = [
            TARGET_VOL,           # vol(m003) + vol(m006)
            TARGET_VOL * 2,       # double
            TARGET_VOL / 2,       # half
            TARGET_VOL * PHI,     # phi multiple
        ]
        
        for t in targets:
            if abs(vol - t) < 0.05:
                candidates.append({
                    "index": i,
                    "name": M.name() if hasattr(M, "name") else f"index_{i}",
                    "volume": vol,
                    "H1": h1,
                    "has_Z5": has_z5,
                    "target": t,
                    "delta": abs(vol - t),
                })
                break
                
    except Exception as e:
        pass
    
    if i % 500 == 0:
        print(f"  Checked {i} manifolds...")

print(f"\nFound {len(candidates)} candidates:")
candidates.sort(key=lambda x: x["delta"])
for c in candidates[:20]:
    z5 = "Z/5" if c["has_Z5"] else "   "
    print(f"  [{c['index']:4d}] vol={c['volume']:.4f}  {z5}  H1={c['H1']:<20}  "
          f"Δ={c['delta']:.5f}  target={c['target']:.4f}")

# Also check exact vol(m003)+vol(m006) candidates with Z/5
print(f"\nZ/5 candidates near {TARGET_VOL:.4f}:")
z5_cands = [c for c in candidates if c["has_Z5"] and abs(c["volume"]-TARGET_VOL)<0.1]
for c in z5_cands:
    print(f"  [{c['index']:4d}] vol={c['volume']:.4f}  H1={c['H1']}")
