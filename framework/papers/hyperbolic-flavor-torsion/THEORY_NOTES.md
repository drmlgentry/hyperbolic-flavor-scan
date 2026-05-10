# Paper 9 Theoretical Development Notes

## Key references found (literature search 2026-03-19)

### Theoretical anchor 1: Ambient Prime Geodesic Theorem
- Balog, Biggs et al., IMRN 2023: "Ambient Prime Geodesic Theorems on Hyperbolic 3-Manifolds"
- https://academic.oup.com/imrn/article/2023/1/588/6381506
- Key result: pi_Gamma(x;J) = |J|/2pi * pi_Gamma(x) + O(e^{5x/3})
  where J is a holonomy interval. This gives equidistribution of twist angles
  as geodesic length grows -- DIRECTLY IMPLIES floors go to zero for ALL classes.
- This is the theoretical proof that phi_floor^(k)(L) -> 0.

### Theoretical anchor 2: Exponential mixing
- Stoyanov (via Dolgopyat): geodesic flow on compact hyperbolic 3-manifolds
  is exponentially mixing. Rate determined by spectral gap lambda_1.
- Li-Pan, Inventiones 2022 (arXiv:2009.12886): exponential mixing for
  geometrically finite case including cusps.
- Implication: mixing rate c relates to lambda_1 via Cheeger-Buser.

### Margulis-Mohammadi-Oh
- Asymptotic for pi_Gamma(x;J) with power savings, for all lattices.
- Dynamical approach. Broader class than arithmetic.

## Proof sketch for Paper 9

**Proposition (from Ambient PGT):**
For each class k in Z/5, phi_floor^(k)(L) -> 0 as L -> infty.

**Proof:** By the Ambient PGT, holonomy is equidistributed as geodesic length grows.
In particular, for any epsilon > 0, the fraction of geodesics with phi_fold < epsilon
in each class k is asymptotically |[0,epsilon]|/pi * (1/5) (equal share per class
since H1 = Z/5 and classes have equal density by abelianization).
Hence every class eventually contains geodesics with phi_fold < epsilon. QED.

**What remains unproven (the conjecture):**
The RATE c_k at which each class achieves phi_floor < epsilon. The Ambient PGT
gives the same asymptotic density for all classes, but the LEADING CONSTANT
controlling the first appearance of near-zero twist elements may differ.
This leading constant is controlled by arithmetic -- specifically the holonomy
representation into PSL(2,C) and the involution k <-> 5-k.

## Remaining development tasks

1. Write Proposition formally using Ambient PGT citation
2. Replace word length with geodesic translation length (use SnapPy length_spectrum)
3. Bootstrap CIs already done (word-length basis)
4. Extend census to m003 for comparison (needs running)
5. Test conjecture: compute trace field of m006 via Snap

## References to add
- Balog-Biggs et al. IMRN 2023
- Margulis (1969, published 2004)
- Stoyanov (2011) for exp mixing on compact hyperbolic 3-manifolds
- Li-Pan Inventiones 2022 arXiv:2009.12886
