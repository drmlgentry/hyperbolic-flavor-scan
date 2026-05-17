"""
tower_comparison.py
Compare M.covers() vs M.low_index_subgroups() for M_PMNS.
This determines whether primes 7 and 29 appear through degree 9.
Run: conda run -n sage python tower_comparison.py
"""
import snappy, re

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
base = {5}

print()
print("METHOD 1: M.covers(n) -- conjugacy classes of subgroups")
seen1 = set()
for n in range(2, 10):
    covers = M.covers(n)
    new = set()
    for C in covers:
        h1 = str(C.homology())
        for m in re.findall(r'Z/(\d+)', h1):
            new |= (full_primes(int(m)) - base - seen1)
        # Also print ALL non-trivial H1
        if h1 not in ('Z/5', '0', 'trivial'):
            orders = re.findall(r'Z/(\d+)', h1)
            if orders:
                primes_here = set()
                for o in orders:
                    primes_here |= full_primes(int(o))
                print(f"  deg={n}: H1={h1}, all_primes={sorted(primes_here)}")
    seen1 |= new
    if new:
        print(f"  *** deg={n}: NEW primes={sorted(new)}")
print(f"\nTotal new primes (covers): {sorted(seen1)}")

print()
print("METHOD 2: M.low_index_subgroups(n) -- ALL subgroups")
seen2 = set()
for n in range(2, 10):
    try:
        subs = M.low_index_subgroups(n)
        new = set()
        all_h1 = []
        for s in subs:
            try:
                h1 = str(s.homology())
                all_h1.append(h1)
                for m in re.findall(r'Z/(\d+)', h1):
                    new |= (full_primes(int(m)) - base - seen2)
            except Exception as e:
                pass
        seen2 |= new
        # Show unique non-trivial H1
        unique = set(h for h in all_h1 if h not in ('Z/5','0','trivial',''))
        for h in sorted(unique):
            orders = re.findall(r'Z/(\d+)', h)
            primes_here = set()
            for o in orders:
                primes_here |= full_primes(int(o))
            print(f"  index={n}: H1={h}, all_primes={sorted(primes_here)}")
        if new:
            print(f"  *** index={n}: NEW primes={sorted(new)}")
        else:
            print(f"  index={n}: no new primes")
    except Exception as e:
        print(f"  index={n}: ERROR {e}")
print(f"\nTotal new primes (low_index): {sorted(seen2)}")

print()
print("CONCLUSION:")
all_primes = sorted(seen1 | seen2)
lucas = {2,3,7,11,29,47}
non_lucas = [p for p in all_primes if p not in lucas]
print(f"  Combined new primes: {all_primes}")
print(f"  Non-Lucas: {non_lucas}")
if 7 in all_primes or 29 in all_primes:
    print("  7 or 29 FOUND -- original claim {2,3,7,11,29} may be correct")
    print("  but depends on which enumeration method is canonical")
else:
    print("  7 and 29 NOT FOUND through degree 9 by either method")
    print("  Original claim {2,3,7,11,29} requires degree > 9")
