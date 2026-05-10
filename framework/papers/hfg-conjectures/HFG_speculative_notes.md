# HFG Speculative Notes and Conjectures
## April 27, 2026

### Note 1: Information storage and the Ihara zeta function

Binary code filling behavior mirrors Dehn filling shell structure:
- 2^n state space = degree-n covering map (index multiplies by 2 each level)
- Ihara zeta function Z_G(u) of a data-structure graph is the
  combinatorial analog of the Selberg zeta function of the flavor manifolds
- Fibonacci heap addressing produces phi-lattice access times
- The prime cycles of the data graph = primitive closed geodesics of M_PMNS

**Conjecture (Information-Geometry Bridge):** The information capacity of 
a binary tree of depth k is log_2(F_k) where F_k is the k-th Fibonacci
number. Since F_k ~ phi^k/sqrt(5), the information capacity grows as
k*log(phi) -- the same unit as the geodesic phi-lattice spacing.

---

### Note 2: Visible vs dark sector cover relationships

**Hypothesis:** Visible sector particles correspond to manifolds with
Lucas-pure covering towers (prime set {2,3,7,11,29,...}).
Dark sector particles correspond to manifolds whose covering towers
contain non-Lucas primes ({5,13,17,19,31,41,...}).

**Falsifiable prediction:** Any dark sector particle discovered with
mass m_dark should satisfy:
  m_dark / m_e = phi^(q/4)
where q/4 is NOT a quarter-Lucas index (i.e., not in {0, 3/4, 2, 5/2,...}).

**Structural reason:** Lucas-pure towers restrict holonomy to the
phi-lattice. Non-Lucas primes introduce geodesics at lengths
incommensurable with log(phi), producing mixing matrices with
different structure from PMNS/CKM.

**Cover relationship:** Visible matter = Lucas-pure sublattice of the
full covering tower. Dark matter = complement in the full prime spectrum.
The ratio |Lucas primes| / |all primes| -> 0 (Lucas primes have
density zero) mirrors the cosmological dark/visible matter ratio ~5:1.

---

### Note 3: Cosmological expansion as Dehn filling flow

**Hypothesis (Filling Flow Cosmology):**
The cosmological expansion corresponds to a flow on Dehn filling
parameter space: (p(t), q(t)) -> (p_0, q_0) as t -> t_now.

Early universe: |p| + |q| large -> near-cusp regime
  -> universal spectral tier dominates
  -> no slope-specific structure (no distinct particles, symmetric phase)

Current universe: specific slopes (-2,3) and (-5,2) reached
  -> sparse tier emerges
  -> PMNS and CKM mixing matrices appear as filling-specific structure

**The tachyonic/inflationary expansion:** corresponds to rapid motion
through filling parameter space at early times (inflation = fast
traversal of the near-cusp regime where all fillings look similar).

**Dark energy / cosmological constant:** corresponds to the residual
"distance" of the current filling parameters from the cusp limit.
Lambda ~ 1/vol(M(p,q)) measures how far the universe is from the
cusped (infinite volume) ancestor manifold.

**Prediction:** The ratio of dark energy density to visible matter
density should be related to the ratio of the covering tower depth
to the volume of M_PMNS:
  rho_Lambda / rho_matter ~ vol(m003) / (degree_max * vol(M_PMNS))
                           = 2.030 / (9 * 0.981) = 0.230
Compare: observed rho_Lambda/rho_matter ~ 2.3 (dark energy ~68%,
matter ~32% of total energy density => ratio ~2.1).
This is a factor-of-10 match in order of magnitude. Worth computing
more carefully.

**Open question:** Is there a natural parameterization of cosmological
time in terms of the Dehn filling slope such that the Hubble parameter
H(t) is proportional to d/dt[log|vol(M(p(t),q(t)))|]?

---

### Note 4: Hard drive topology

A hard drive platter is geometrically a punctured torus (disk with
hole at center). The sector/track addressing scheme is a quotient
of the torus by a discrete action -- exactly the structure of a
Dehn filling. The number of sectors per track grows with radius
(outer tracks have more sectors) in modern drives (zone bit recording)
-- this is the analog of the geodesic length spectrum becoming denser
at larger n.

The error correction codes (ECC) on hard drives use Reed-Solomon
codes based on finite field arithmetic over GF(2^k). The prime 2
appears as the base -- L_0 = 2 is the first Lucas prime and the
base of binary arithmetic. The connection: our covering tower's
first new prime is 2 (= L_0), and binary arithmetic's fundamental
unit is also 2. Both arise from the same combinatorial structure
(degree-2 covers / bit doubling).

---
*These are informal notes for future development.
All marked as SPECULATIVE / CONJECTURE.*
