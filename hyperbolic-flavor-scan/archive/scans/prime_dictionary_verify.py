import snappy
import re

M3 = snappy.OrientableClosedCensus[1]
v3 = float(M3.volume())

# Extended predicted prime set with geometric justification
predicted = {2, 3, 5, 7, 11, 13, 29}
print("Extended predicted prime set:", sorted(predicted))
print("Sources:")
print("  2  = |p_PMNS|")
print("  3  = |q_PMNS|")
print("  5  = H1(m003) torsion")
print("  7  = |p_PMNS| + |p_CKM|")
print("  11 = det(PMNS slope, CKM slope)")
print("  13 = ||PMNS||^2 = 4+9")
print("  29 = ||CKM||^2  = 25+4")

def prime_factors(n):
    n = abs(n)
    if n <= 1:
        return set()
    factors = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    return factors

all_unexpected = {}
total_covers = 0
v3_exact = v3

print(f"\nScanning degrees 2-6 (tolerance 1e-5 for cleaner matching)...")
for deg in range(2, 7):
    target = deg * v3_exact
    covers = []
    for i in range(len(snappy.OrientableClosedCensus)):
        M = snappy.OrientableClosedCensus[i]
        v = float(M.volume())
        if abs(v - target) < 1e-5:   # tighter tolerance
            covers.append((i, M.name(), v, str(M.homology())))
    
    print(f"\n--- Degree {deg} ({len(covers)} covers, vol={target:.8f}) ---")
    for idx, name, vol, h1 in covers:
        nums = [int(x) for x in re.findall(r'\d+', h1) if int(x) > 1]
        primes = set()
        for num in nums:
            primes |= prime_factors(num)
        unexpected = primes - predicted
        if unexpected:
            all_unexpected[idx] = (name, h1, unexpected)
            print(f"  census[{idx:5d}] {name:10s} H1={h1:30s} UNEXPECTED: {unexpected}")
        else:
            print(f"  census[{idx:5d}] {name:10s} H1={h1:30s} ok {sorted(primes)}")
        total_covers += 1

print(f"\n=== Summary ===")
print(f"Total covers scanned: {total_covers}")
print(f"Manifolds with primes outside predicted set: {len(all_unexpected)}")
if all_unexpected:
    for idx, (name, h1, unexp) in all_unexpected.items():
        print(f"  census[{idx}] {name} H1={h1} unexpected primes: {unexp}")
    print("\nThese require explanation or indicate prediction needs refinement.")
else:
    print("ALL primes in {2,3,5,7,11,13,29} -- conjecture holds through degree 6")
