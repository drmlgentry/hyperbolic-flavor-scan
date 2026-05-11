## Priority 3: CP Phase from A-Factor — RESULT: NOT FOUND (May 10, 2026)

### Methods tested
1. Twist angles of individual word representatives -- UNSTABLE (std 36 deg)
   Best hit: AAAb/aaBa/aaaB/bAAA all give phi=200.649 deg in default triangulation
   but std=35-37 deg across 20 random triangulations -- triangulation artifact
   
2. Manifold-level topological invariants -- STABLE but wrong values
   - Dedekind sum D(-2,3) * 360 = 20.0 deg
   - Chern-Simons * 360 = 90.0 deg (CS=1/4 confirmed separately)
   - Rademacher Phi(-2,3) = 0.667 -> 240 deg
   - Slope complex phase arg(-2+3i) = 123.69 deg

### Interesting observations
- Filling slope ratio arctan(2/3) = 33.69 deg ~ theta_12 = 33.65 deg (diff 0.04 deg)
  This may be worth noting but is not the CP phase
- The CP phase delta_CP = 197 deg does not match any triangulation-independent 
  quantity we can compute from m003(-2,3)

### Conclusion
Priority 3 is OPEN. The CP phase is not currently derivable from m003(-2,3)
by any triangulation-independent geometric method. The J_PMNS != 0 result
from the homology no-go theorem remains valid (it proves CP violation exists
in the lepton sector) but does not give the specific phase value.

Possible future approaches:
- Rho invariant of the rational homology sphere (requires more topology)
- Secondary invariants of the flat PSL(2,C) bundle over m003
- Relationship between delta_CP and the character variety coordinates
