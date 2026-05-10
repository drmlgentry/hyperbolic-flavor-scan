import numpy as np

# Filling coefficients
pmns = np.array([-2, 3])
ckm  = np.array([-5, 2])

# Intersection number on boundary torus = |det([pmns, ckm])|
det = int(pmns[0]*ckm[1] - pmns[1]*ckm[0])
print(f"PMNS filling slope: {pmns}")
print(f"CKM  filling slope: {ckm}")
print(f"Intersection number |det|: {abs(det)}")
print(f"  = |({pmns[0]})({ckm[1]}) - ({pmns[1]})({ckm[0]})| = |{pmns[0]*ckm[1]} - {pmns[1]*ckm[0]}|")

# Norms
pmns_norm = np.sqrt(pmns @ pmns)
ckm_norm  = np.sqrt(ckm @ ckm)
print(f"\n||PMNS slope|| = sqrt({pmns@pmns}) = {pmns_norm:.6f}")
print(f"||CKM  slope|| = sqrt({ckm@ckm})  = {ckm_norm:.6f}")
print(f"Ratio of norms: {ckm_norm/pmns_norm:.6f}")

# Angle between slopes
cos_theta = (pmns @ ckm) / (pmns_norm * ckm_norm)
theta = np.arccos(np.clip(cos_theta, -1, 1))
print(f"\nAngle between slopes: {np.degrees(theta):.4f} degrees")
print(f"                      {theta:.6f} radians")

# Farey/continued fraction analysis of slopes
from fractions import Fraction
pmns_slope = Fraction(pmns[0], pmns[1])
ckm_slope  = Fraction(ckm[0],  ckm[1])
print(f"\nPMNS slope as fraction: {pmns_slope} = {float(pmns_slope):.6f}")
print(f"CKM  slope as fraction: {ckm_slope}  = {float(ckm_slope):.6f}")
print(f"Sum of slopes: {pmns_slope + ckm_slope} = {float(pmns_slope + ckm_slope):.6f}")
print(f"Product of slopes: {pmns_slope * ckm_slope} = {float(pmns_slope * ckm_slope):.6f}")

# Farey distance
farey_dist = abs(pmns[0]*ckm[1] - pmns[1]*ckm[0])
print(f"\nFarey distance (= intersection number): {farey_dist}")
print(f"Primes in play: factorize {farey_dist} = ", end="")
n = farey_dist
factors = []
d = 2
while d * d <= n:
    while n % d == 0:
        factors.append(d)
        n //= d
    d += 1
if n > 1:
    factors.append(n)
print(" × ".join(map(str, factors)) if factors else "1")

# Check: does 7 appear from any natural combination?
print(f"\nChecking for 7 in filling geometry:")
print(f"  |pmns| + |ckm| components: {list(abs(pmns))} + {list(abs(ckm))}")
print(f"  Sum of |p| values: {abs(pmns[0]) + abs(ckm[0])} = {abs(pmns[0])}+{abs(ckm[0])}")
print(f"  Sum of |q| values: {abs(pmns[1]) + abs(ckm[1])} = {abs(pmns[1])}+{abs(ckm[1])}")
print(f"  |p_PMNS| + |q_CKM|: {abs(pmns[0]) + abs(ckm[1])}")
print(f"  |q_PMNS| + |p_CKM|: {abs(pmns[1]) + abs(ckm[0])}")
