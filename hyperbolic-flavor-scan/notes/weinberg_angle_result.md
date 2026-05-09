# Weinberg Angle Result — May 8, 2026

## Discovery
Pairwise angle between AAb and BaaB on m006 = 28.1790 deg
Weinberg angle theta_W = 28.1826 deg
Deviation: 0.0036 deg (0.013%)

## BaaB is isospectral with the CKM triple (confirmed)
|tr(aaB)| = |tr(AbA)| = |tr(AAb)| = |tr(BaaB)| = |tr(Baa)| = 2.014574
tr(BaaB) = conj(tr(AAb)) -- orientation reversal of m006 (amphicheiral)

## Interpretation
BaaB = b^{-1} * a^2 * b^{-1} is the orientation-reversed conjugate of AAb.
m006 is amphicheiral: admits orientation-reversing isometry.
The angle between AAb and its mirror BaaB = theta_W to 0.013%.

## Physical meaning (conjecture)
Weinberg angle = geometric angle between a CKM holonomy word
and its orientation-reversed mirror in the amphicheiral CKM manifold.
Connects electroweak mixing angle to amphicheirality of m006.

## Also found in m003 (PMNS)
BAB vs BBAb: angle = 28.1937 deg (diff 0.011 deg from theta_W)

## Status: CONJECTURE
Needs derivation from first principles showing angle = theta_W exactly.

## Next steps
1. Verify axis angle between AAb and BaaB directly (confirm sweep result)
2. Check Weinberg angle in degree-2 cover (index 39) and degree-3 (index 238)
3. Check if angle is invariant across Dehn fillings

## Dehn filling sweep (May 8, 2026)
AAb vs BaaB angle across Dehn fillings of cusped m006:
  (-5,2)  vol=2.028853  angle=28.1790  diff=-0.0036  *** CKM filling ***
  (-4,2)  vol=1.843586  angle= 8.2663  diff=-19.92
  (-6,2)  vol=2.173118  angle=15.9485  diff=-12.23
  (-5,3)  vol=2.184756  angle=20.7897  diff= -7.39
  (-5,1)  vol=1.941503  angle=16.0174  diff=-12.17

Covering tower of m003:
  idx39  (deg-2): AAb/BaaB=60.25, BAB/BBAb=39.33 -- NO match
  idx238 (deg-3): AAb/BaaB=54.88, BAB/BBAb=11.79 -- NO match

## CONCLUSION
theta_W is NOT a universal property of the cusped ancestor.
It is specific to the EXACT Dehn fillings (-5,2) and (-2,3) that
optimize CKM and PMNS respectively.
The weak mixing angle is a geometric selection -- encoded precisely
in the filling slopes, not in the topology alone.

## Status: RESULT (numerical, needs derivation)

## Axis geometry — May 8, 2026 (second computation)

### All Weinberg-angle axes lie in the xz-plane (y=0 exactly)
Reflects a Z/2 symmetry of the holonomy representation of both m006 and m003.

### The dot product identity
|n(AAb) . n(BaaB)| = 0.881521
cos(theta_W) = M_W/M_Z = 0.881447
Difference: 0.000075  (agreement to 4 sig figs in the cosine)

### The separation angle in both manifolds
m006: separation(AAb, BaaB) in xz-plane = 208.1815 deg
m003: separation(BAB, BBAb) in xz-plane = 208.1914 deg
Difference between manifolds: 0.0099 deg

### The compact statement
|n_AAb . n_BaaB| = cos(theta_W) = M_W/M_Z
in the invariant xz-plane of the holonomy of m006(-5,2).

Same identity holds in m003(-2,3) with words BAB/BBAb.

### Status: NUMERICAL THEOREM
The identity |dot product| = M_W/M_Z holds to 4 significant figures.
Needs algebraic derivation from the arithmetic of the filling slopes.

## RETRACTION OF WEINBERG ANGLE CLAIM — May 8, 2026

The dot product |n(AAb) . n(BaaB)| = cos(theta_W) is TRIANGULATION-DEPENDENT.

Across 20 random triangulations of m006(-5,2), the dot product ranges
from 0.054 to 0.994 (std=0.29). The value 0.881 is not distinguished.

The result was an artifact of SnapPy's default triangulation, not a
geometric property of the manifold.

STATUS: Demoted from result to artifact. Do not include in papers.

What remains valid:
- The xz-plane invariance (y=0) is also likely triangulation-dependent
  and should be checked independently
- The Dehn surgery ranking (#1 out of 87) is also likely affected
  since it uses the same triangulation-dependent axis extraction
- The CKM/PMNS fitness results (0.017, 0.019) ARE valid because
  the fitness metric compares mixing angles, not axis dot products --
  these are representation-independent observables

What to investigate instead:
- Find a triangulation-INDEPENDENT quantity that encodes theta_W
- The character variety chi(pi_1(m006)) -> C is intrinsic
- Trace functions tr(rho(g)) are conjugacy-invariant
- Could theta_W appear as a ratio of trace lengths or eigenvalues?

## SECOND RETRACTION: xz-plane invariance — May 8, 2026

The y=0 result was a mathematical tautology, not a geometric theorem.

BUG: In axis_imlogm_v1, ny is computed as:
  imL = L.imag          # imL is a real-valued matrix
  ny = (1j*(imL[1,0] - imL[0,1])).real

Since imL[1,0] and imL[0,1] are real numbers,
  1j * (real - real) = purely imaginary number
  .real = 0  ALWAYS

This is a mathematical identity -- ny=0 for ANY matrix, ANY manifold.
The "xz-plane invariance theorem" was a tautology built into the formula.

The correct ny formula is:
  ny = (1j*(L[1,0] - L[0,1])).real   # L is the full complex logm
This gives non-zero y-components as confirmed by the second computation.

## SUMMARY OF RETRACTIONS TODAY

1. Weinberg angle result: triangulation-dependent artifact
2. xz-plane invariance: mathematical tautology in axis formula

## WHAT REMAINS VALID

All results that do NOT use the imlogm axis formula:
- CKM/PMNS fitness (0.017, 0.019) -- uses pairwise angles between axes
  but the FITNESS VALUE is computed from mixing matrix entries, not
  raw axis dot products. Need to check if fitness is also affected.
- Covering tower, Lucas structure, CS invariants -- do not use axes
- Volume ratios, partition functions -- do not use axes

## IMMEDIATE ACTION NEEDED

Verify that the CKM and PMNS FITNESS COMPUTATIONS are not affected
by the v1 formula bug. The fitness uses QR decomposition of the
overlap matrix, which depends on the axis vectors. If the axis
vectors have wrong y-components due to the bug, the fitness values
may be artifacts too.
