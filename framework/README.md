# Hyperbolic Flavor Geometry

A research program deriving Standard Model flavor parameters from the geometry
of compact hyperbolic 3-manifolds via the Iwasawa decomposition PSL(2,C) = KAN.

## Core Result

The Iwasawa decomposition distributes SM flavor structure across two manifolds:
- **K-factor** → CKM quark mixing | manifold **m006** (vol=2.029, H1=Z/5) | fitness F=0.01729
- **N-factor** → PMNS lepton mixing | manifold **m003** (Meyerhoff, vol=0.981, H1=Z/5) | fitness F=0.01897
- **A-factor** → CP phases and full flavor spectrum from loxodromic twist angles φ(γ)

Both optimal manifolds have H1=Z/5 (odd-torsion conjecture).
m003 is the **Meyerhoff manifold** (OrientableClosedCensus[1]); distinct from the
Weeks manifold (OrientableClosedCensus[0], vol=0.9427).

## Submission Portfolio (8 papers)

| # | Title | Journal | Status | ID |
|---|---|---|---|---|
| 1 | CP Phases from Holonomy | Nuclear Physics B | Under review | NPB-S-26-00539 |
| 2 | Hyperbolic Log Lattices | Geometriae Dedicata | Under review | rs-9071491 |
| 3 | Discrete Mixing Operators | Nuclear Physics B | Under review | NPB-S-26-00540 |
| 4 | Scale-Free Quadratic Forms | Nuclear Physics B | Under review | NPB-S-26-00538 |
| 5 | CKM from Geodesic Axes | Physical Review D | Under review | es2026mar11_966 |
| 6 | PMNS from Borel Structure | Physical Review D | Under review | es2026mar13_942 |
| 7 | CP Violation from A-Factor | Physical Review D | Under review | es2026mar14_722 |
| 8 | Twist Angle Spectrum | JHEP | Resubmitting (NPB desk-rejected: out of scope) | — |

## Key Numerical Results

### CKM (m006)
- Words: aaB/AbA/AAb | fitness=0.01729 | J=0
- Axis angles: θ12=48.2°, θ13=77.5°, θ23=68.4°
- Triple isospectrality: φ≈92.49° for all three words → explains J_CKM≈0

### PMNS (m003 = Meyerhoff manifold)
- Words: aa/ab/aB | column perm (0,2,1) | fitness=0.01897
- CP phase: φ_aa − φ_ab + φ_aB = 203.5° vs PDG 197° (3.3% error, 0 free params)

### Twist Angle Spectrum (Paper 8)
- m006 word aa: 180°−φ = 67.65° ~ δ_CKM=68.0° (0.51% error)
- m006 word abbAB: 180°−φ = 33.62° ~ θ12_ν=33.41° (0.63% error)
- m003 ratio |λ(bbbb)|/|λ(bAbA)| = 3.2910 ~ m̄_b/m̄_c = 3.2913 (0.01%, MS-bar)
- Length-7 census: spectral gap [0.005°, 1.61°] on m006; θ13_CKM sits inside gap
- Selberg estimates: λ1(m003)≈2.48, λ1(m006)≈2.82 (first for these manifolds)

## Repository Structure
```
papers/
  hyperbolic-flavor-ckm/      # Paper 5 (PRD) - CKM from geodesic axes
  hyperbolic-flavor-pmns/     # Paper 6 (PRD) - PMNS from Borel structure
  hyperbolic-flavor-cp/       # Paper 7 (PRD) - CP violation from A-factor
  hyperbolic-flavor-twist/    # Paper 8 (JHEP) - Twist angle spectrum
code/                         # Utility scripts
```

Scan scripts and census data: https://github.com/drmlgentry/hyperbolic-flavor-scan

## Author

Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663
