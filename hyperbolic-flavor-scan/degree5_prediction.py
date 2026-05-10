import snappy
import numpy as np

# Established prime dictionary
# deg 1: {5}
# deg 2: {5, 7, 11}
# deg 3: {2, 3, 5, 7}
# deg 4: {2, 3, 5, 7}
# Prediction for deg 5: no new primes beyond {2,3,5,7,11}
# Falsification criterion: any prime > 11 appearing in deg-5 cover H1

M3 = snappy.OrientableClosedCensus[1]
v3 = float(M3.volume())
target = 5 * v3

print(f"M3/PMNS vol={v3:.8f}")
print(f"Degree-5 target vol={target:.8f}")
print(f"Searching full census for degree-5 covers...\n")

deg5 = []
n = len(snappy.OrientableClosedCensus)
for i in range(n):
    if i % 2000 == 0:
        print(f"  Progress: {i}/{n}...")
    M = snappy.OrientableClosedCensus[i]
    v = float(M.volume())
    if abs(v - target) < 1e-4:
        deg5.append((i, M.name(), v, str(M.homology())))

print(f"\n=== Degree-5 covers of M3/PMNS (vol={target:.6f}) ===")
print(f"Count: {len(deg5)}")
for idx, name, vol, h1 in deg5:
    print(f"  census[{idx:5d}] {name:10s} vol={vol:.8f} H1={h1}")

# Extract all primes appearing and check against prediction
def prime_factors(n):
    n = abs(n)
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

import re
predicted = {2, 3, 5, 7, 11}
new_primes_found = set()

print(f"\n=== Prime Analysis ===")
print(f"Predicted prime set: {sorted(predicted)}")
for idx, name, vol, h1 in deg5:
    # Parse H1 string for integers
    nums = [int(x) for x in re.findall(r'\d+', h1)]
    primes = set()
    for num in nums:
        primes |= prime_factors(num)
    unexpected = primes - predicted
    new_primes_found |= unexpected
    status = "UNEXPECTED: " + str(unexpected) if unexpected else "ok"
    print(f"  census[{idx:5d}] H1={h1:25s} primes={sorted(primes)} {status}")

print(f"\nNew primes beyond prediction: {sorted(new_primes_found) if new_primes_found else 'NONE'}")
if not new_primes_found:
    print("PREDICTION CONFIRMED: degree-5 tower contains only primes in {2,3,5,7,11}")
else:
    print("PREDICTION FALSIFIED: new primes found, pattern breaks at degree 5")

# Also: natural prediction for degree-5 from slope arithmetic
pmns = np.array([-2, 3])
ckm  = np.array([-5, 2])
print(f"\n=== What degree-5 arithmetic would predict ===")
print(f"5*|p_PMNS| = {5*abs(pmns[0])} = 2*5 (not prime)")
print(f"5*|q_PMNS| = {5*abs(pmns[1])} = 3*5 (not prime)")
print(f"|p|^2 + |q|^2 = {pmns[0]**2 + pmns[1]**2} <- 13 (prime, first appeared degree-2 norm)")
print(f"det + |p| + |q| = {11 + 2 + 3} = 16 (not prime)")
print(f"det - |p|*|q| = {11 - 2*3} = 5 (known)")
print(f"Prediction: no prime beyond 11 should appear at degree 5")
