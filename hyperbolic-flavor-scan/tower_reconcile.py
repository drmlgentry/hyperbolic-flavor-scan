"""
tower_reconcile.py
==================
Compare two methods of finding covers of M_PMNS:
  Method A: M.covers(degree) -- algebraic, finds ALL covers
  Method B: volume matching in OrientableClosedCensus -- finds census covers

Question: do these find the same covers?
If not: the Dehn filling slopes paper's prime dictionary
may describe a different structure than the Lucas-pure tower.

Run: conda run -n sage python tower_reconcile.py
"""
import snappy
import numpy as np

M_PMNS = snappy.OrientableClosedCensus[1]
vol_PMNS = float(M_PMNS.volume())

print("="*65)
print("COVERING TOWER RECONCILIATION")
print("="*65)
print(f"\nM_PMNS = {M_PMNS.name()}, vol={vol_PMNS:.8f}, H1={M_PMNS.homology()}")
print()

for degree in range(2, 7):
    target_vol = degree * vol_PMNS
    tol = 1e-4

    print(f"{'='*55}")
    print(f"DEGREE {degree}  (target vol = {target_vol:.6f})")
    print(f"{'='*55}")

    # Method A: algebraic
    print(f"\nMethod A: M.covers({degree})")
    try:
        covers_A = M_PMNS.covers(degree)
        print(f"  Found {len(covers_A)} covers")
        for i,C in enumerate(covers_A):
            try:
                h1 = C.homology()
                vol = float(C.volume())
                print(f"  [{i}] vol={vol:.6f}, H1={h1}")
            except Exception as e:
                print(f"  [{i}] error: {e}")
    except Exception as e:
        print(f"  Error: {e}")

    # Method B: census volume matching
    print(f"\nMethod B: Census volume matching (tol={tol})")
    census_hits = []
    for idx in range(0, 500):
        try:
            M_c = snappy.OrientableClosedCensus[idx]
            vol_c = float(M_c.volume())
            if abs(vol_c - target_vol) < tol:
                h1 = M_c.homology()
                census_hits.append((idx, M_c.name(), vol_c, h1))
        except:
            break

    if census_hits:
        print(f"  Found {len(census_hits)} in census:")
        for idx, name, vol, h1 in census_hits:
            print(f"  [idx={idx:>4}] {name}, vol={vol:.6f}, H1={h1}")
    else:
        print(f"  None found in first 500 census entries")

    # Compare H1 prime sets
    print()

print("="*65)
print("PRIME ANALYSIS: what primes appear in H1?")
print("="*65)
print()

all_primes_A = set()
all_primes_B = set()

def get_h1_primes(h1_str):
    """Extract prime factors from H1 string."""
    import re
    primes = set()
    # Find all numbers in the H1 description
    nums = re.findall(r'\d+', str(h1_str))
    for n in nums:
        n = int(n)
        if n <= 1: continue
        # Factor n
        for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]:
            while n % p == 0:
                primes.add(p)
                n //= p
        if n > 1:
            primes.add(n)
    return primes

print("Algebraic covers (Method A):")
for degree in range(2, 7):
    try:
        covers = M_PMNS.covers(degree)
        for C in covers:
            try:
                h1 = C.homology()
                p = get_h1_primes(h1)
                all_primes_A |= p
                print(f"  deg={degree}: H1={h1}, primes={sorted(p)}")
            except:
                pass
    except:
        pass

print(f"\nAll new primes in algebraic tower (deg 2-6): {sorted(all_primes_A)}")

print("\nCensus covers (Method B):")
for degree in range(2, 7):
    target_vol = degree * vol_PMNS
    for idx in range(0, 500):
        try:
            M_c = snappy.OrientableClosedCensus[idx]
            vol_c = float(M_c.volume())
            if abs(vol_c - target_vol) < 1e-4:
                h1 = M_c.homology()
                p = get_h1_primes(h1)
                all_primes_B |= p
                print(f"  deg={degree}, idx={idx}: H1={h1}, primes={sorted(p)}")
        except:
            break

# Remove primes from M_PMNS itself
base_primes = get_h1_primes(M_PMNS.homology())
new_primes_A = all_primes_A - base_primes
new_primes_B = all_primes_B - base_primes

print(f"\nBase H1 primes of M_PMNS: {sorted(base_primes)}")
print(f"\nNEW primes (algebraic covers): {sorted(new_primes_A)}")
print(f"NEW primes (census covers):    {sorted(new_primes_B)}")
print()
print("Are they the same?", sorted(new_primes_A) == sorted(new_primes_B))
print()
if new_primes_A != new_primes_B:
    print("In algebraic but not census:", sorted(new_primes_A - new_primes_B))
    print("In census but not algebraic:", sorted(new_primes_B - new_primes_A))
    print()
    print("INTERPRETATION:")
    print("  If census has extra primes: census includes non-covering")
    print("  manifolds of same volume (coincidental, not covers)")
    print("  If algebraic has extra primes: covers exist outside census")
print()
print("="*65)
print("LUCAS CHECK")
print("="*65)
print()
lucas = [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322, 521]
lucas_primes = {x for x in lucas if x > 1 and all(x%i!=0 for i in range(2,x)) and x > 1}
print(f"Lucas prime numbers: {sorted(lucas_primes)}")
print(f"New algebraic primes: {sorted(new_primes_A)}")
print(f"Algebraic primes are Lucas primes: "
      f"{new_primes_A.issubset(lucas_primes)}")
print(f"New census primes: {sorted(new_primes_B)}")
print(f"Census primes are Lucas primes: "
      f"{new_primes_B.issubset(lucas_primes)}")
