# Compute Chern-Simons and eta invariants for HFG manifolds
# Run in WSL conda sage environment

import snappy
from snappy import Manifold
import cypari2

pari = cypari2.Pari()

manifolds = {
    'M_PMNS = m003(-2,3)': 'm003(-2,3)',
    'M_CKM  = m006(-5,2)': 'm006(-5,2)',
    'm019 (SU(4) parent)': 'm019',
    'm178':               'm178',
    'm179':               'm179',
    'v1024':              'v1024',
    't03293':             't03293',
    'v2603':              'v2603',
}

v0 = Manifold('m003(-2,3)').volume()
print(f"v0 = {float(v0):.10f}")
print()
print(f"{'Manifold':<30} {'Volume':>14} {'vol/v0':>10} {'CS (mod 1)':>14} {'CS*v0':>12}")
print("-"*85)

for name, mfld in manifolds.items():
    M = Manifold(mfld)
    vol = float(M.volume())
    ratio = vol / float(v0)
    try:
        cs = float(M.chern_simons())
        cs_v0 = cs * float(v0)
    except:
        cs = None
        cs_v0 = None
    
    cs_str = f"{cs:.8f}" if cs is not None else "N/A"
    csv0_str = f"{cs_v0:.6f}" if cs_v0 is not None else "N/A"
    print(f"{name:<30} {vol:>14.10f} {ratio:>10.6f} {cs_str:>14} {csv0_str:>12}")

print()
print("# Checking CS relationships")
print("# If CS(M) = p/q * v0 for rational p/q, that's a new invariant")
print()

# Also check cusped parents
print("Cusped parents:")
for name, mfld in [('m003','m003'),('m006','m006'),('m019','m019')]:
    M = Manifold(mfld)
    try:
        cs = M.chern_simons()
        print(f"  {name}: CS = {cs} (cusped -- may be 0 or undefined)")
    except Exception as e:
        print(f"  {name}: CS error: {e}")

