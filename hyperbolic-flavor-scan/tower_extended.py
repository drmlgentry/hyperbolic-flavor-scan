"""
tower_extended.py
Extend covering tower computation beyond degree 9 to find primes 7 and 29.
Optimized: early exit once 7 and 29 found, progress reporting.
Run: conda run -n sage python tower_extended.py
"""
import snappy, re, time

def full_primes(n):
    primes = set()
    temp = n
    p = 2
    while p*p <= temp:
        if temp % p == 0:
            primes.add(p)
            while temp % p == 0:
                temp //= p
        p += 1
    if temp > 1:
        primes.add(temp)
    return primes

M = snappy.OrientableClosedCensus[1]
print(f"M_PMNS: {M.name()}, H1={M.homology()}")
print("Extending covering tower beyond degree 9...")
print("Target: find degrees where primes 7 and 29 first appear.")
print()

base = {5}
seen = set()
targets = {7, 29}
found = {}

MAX_DEG = 20  # extend well beyond 9

for deg in range(2, MAX_DEG + 1):
    t0 = time.time()
    try:
        covers = M.covers(deg)
        n_covers = len(covers)
    except Exception as e:
        print(f"deg={deg}: ERROR {e}")
        continue

    new_at_deg = set()
    for C in covers:
        h1 = str(C.homology())
        for m in re.findall(r'Z/(\d+)', h1):
            new_at_deg |= (full_primes(int(m)) - base - seen)

    seen |= new_at_deg
    elapsed = time.time() - t0

    # Report
    if new_at_deg:
        lucas = {2,3,7,11,29,47,199,3571}
        non_luc = [p for p in new_at_deg if p not in lucas]
        print(f"deg={deg:2d}: {n_covers:4d} covers, "
              f"NEW={sorted(new_at_deg)}, non-Lucas={non_luc}, "
              f"t={elapsed:.1f}s")
        # Record first appearance
        for p in new_at_deg:
            if p in targets and p not in found:
                found[p] = deg
                print(f"  *** PRIME {p} = L_{[2,3,7,11,29,47].index(p)*2 if p in [2,3,7,11,29,47] else '?'} FOUND at degree {deg} ***")
    else:
        print(f"deg={deg:2d}: {n_covers:4d} covers, (none), t={elapsed:.1f}s")

    # Early exit if both targets found
    if targets <= set(found.keys()):
        print(f"\nBoth targets found! Stopping at degree {deg}.")
        break

    # Also stop if taking too long per degree
    if elapsed > 120:
        print(f"\nStopping: degree {deg} took {elapsed:.0f}s (too slow for higher degrees)")
        break

print()
print("="*55)
print("SUMMARY")
print("="*55)
print(f"All new primes found: {sorted(seen)}")
lucas_names = {2:'L_0', 3:'L_2', 7:'L_4', 11:'L_5', 29:'L_7', 47:'L_8'}
for p in sorted(seen):
    name = lucas_names.get(p, f"p={p}")
    deg = found.get(p, '?')
    print(f"  {p} ({name}): first appears at degree {deg}")
print()
if 7 not in seen:
    print("Prime 7 NOT found through computed degrees.")
    print("Original claim requires computation beyond current degree.")
if 29 not in seen:
    print("Prime 29 NOT found through computed degrees.")
non_lucas = [p for p in seen if p not in lucas_names]
if non_lucas:
    print(f"NON-LUCAS PRIMES: {non_lucas} -- LUCAS PURITY VIOLATED")
else:
    print("Lucas-purity holds for all primes found.")
