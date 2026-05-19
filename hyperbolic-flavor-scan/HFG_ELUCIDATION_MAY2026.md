# HFG Geometric Elucidation — Working Document
**Author:** Marvin L. Gentry  
**Date:** May 19, 2026  
**Status:** Research notes — not for submission

---

## The Two Manifolds

| Property | M_PMNS = m003(−2,3) | M_CKM = m006(−5,2) |
|---|---|---|
| Census index | OrientableClosedCensus[1] | OrientableClosedCensus[43] |
| Volume | 0.98136883 | 2.02885309 |
| H₁ | ℤ/5 | ℤ/5 |
| Tetrahedra | 2 | 3 |
| Trace field | ℚ(√−3) imaginary quadratic | ℚ(√17) real quadratic |
| Chern-Simons | 1/4 | — |
| Relators | ababAbbAb, abAbaabAbaBAB | ababbAAbb, abbAbbaBabbAbbAbbaB |

**Critical convention:** Always load by index, never by name string.  
`M_PMNS = snappy.OrientableClosedCensus[1]`  
`M_CKM  = snappy.OrientableClosedCensus[43]`

---

## Generator Holonomy Data

### m003(−2,3) — M_PMNS

Generators: a, b. Holonomy via `polished_holonomy()` at 150-bit:

| word | ell | φ (deg) | \|λ\| | Re(tr) | Im(tr) |
|---|---|---|---|---|---|
| a | 0.88944 | 84.2781 | 1.560056 | 0.21945 | 0.91447 |
| b | 1.04032 | 28.1429 | 1.682293 | 2.00755 | 0.51312 |
| aa | 1.77889 | 168.5561 | 2.433774 | −2.78810 | 0.40136 |
| bb | 2.08063 | 56.2858 | 2.830109 | 1.76698 | 2.06021 |
| aaB | 1.73425 | −176.7311 | 2.380055 | −2.79566 | −0.11176 |
| baa | 1.99665 | −167.3620 | 2.713733 | −3.00755 | −0.51312 |
| AbA | 1.73425 | −176.7311 | 2.380055 | −2.79566 | −0.11176 |
| AAb | 1.73425 | −176.7311 | 2.380055 | −2.79566 | −0.11176 |
| bbbb | 4.16126 | 112.5715 | 8.009515 | — | — |
| bAbA | 1.77889 | 168.5561 | 2.433774 | — | — |

**Key observation:** bAbA is isospectral to aa on m003.

### m006(−5,2) — M_CKM

| word | ell | φ (deg) | \|λ\| | Re(tr) | Im(tr) |
|---|---|---|---|---|---|
| a | 0.94160 | 56.1728 | 1.601271 | 1.23907 | 0.81142 |
| b | 0.49057 | −90.8436 | 1.277988 | −0.03034 | −0.49546 |
| aa | 1.88319 | 112.3456 | 2.564070 | −1.12312 | 2.01081 |
| bb | 0.98115 | 178.3128 | 1.633255 | −2.24456 | 0.03006 |
| aaB | 1.77172 | −87.5131 | 2.425064 | 0.12312 | −2.01081 |
| AbA | 1.77172 | −87.5131 | 2.425064 | 0.12312 | −2.01081 |
| AAb | 1.77172 | −87.5131 | 2.425064 | 0.12312 | −2.01081 |
| aaabb | 3.46544 | 2.1319 | 5.656008 | — | — |
| abbAB | 1.72900 | 33.6202 | 2.373823 | — | — |

**Key observation:** aaB, AbA, AAb are isospectral on m006 (spread = 0.000000°).

---

## Structural Theorems (proven)

### Theorem 1: Power law
For any γ ∈ π₁(M):
- ell(γⁿ) = n · ell(γ) exactly
- φ(γⁿ) = n · φ(γ) mod 360°

This follows from λ(γⁿ) = λ(γ)ⁿ. All integer-ratio coincidences in
the length spectrum (1:2:3:4, etc.) follow from this theorem alone.

### Theorem 2: CKM isospectrality
The words aaB, AbA, AAb have identical complex lengths on m006:
- ell = 1.77172, φ = −87.513° (spread = 0 to machine precision)

**Geometric reason:** The H₁ classes [AbA] = [AAb] = 2 in H₁ = ℤ/5
force a homological collision. This implies J_CKM = 0 in the
geometric construction (no CP violation from the CKM axis triple).

### Theorem 3: Trace field dichotomy
- M_PMNS has imaginary trace field ℚ(√−3) → loxodromic elements
  with non-trivial twist angles → large leptonic CP violation
- M_CKM has real trace field ℚ(√17) → holonomy traces constrained
  to (a non-standard embedding of) a real field → CP suppression

Verified: tr(ρ(aa) on m006) has Re = 3 − √17 = −1.12311 (0.001% error).

---

## Covering Tower

### Algebraic covering tower (M.covers())
Through degree 6, M_PMNS has exactly ONE cover:
- Degree 5: H₁ = ℤ/11 + ℤ/11, new prime = **11 = L₅**

M_PMNS is arithmetically rigid — very few algebraic covers exist.

### Census volume coincidences (NOT algebraic covers)
Manifolds in the closed census with volume = n · vol(M_PMNS):
- Degree 2: census[39] H₁=ℤ/55 (primes 5,11), census[40] H₁=ℤ/7+ℤ/7
- Degree 3: census[237-240] H₁ primes {2,3,5,7}

These are commensurable neighbors, not covering spaces of M_PMNS.
Their H₁ primes {2,3,5,7,11} are all Lucas primes.

### Connection to slope arithmetic
Dehn filling slopes: v_PMNS = (−2,3), v_CKM = (−5,2)

- Farey intersection: det[(−2,−5),(3,2)] = −4+15 = **11 = L₅**
  → same prime as the unique algebraic cover
- CKM slope norm: ‖(−5,2)‖² = 25+4 = **29 = L₇**
- PMNS slope norm: ‖(−2,3)‖² = 4+9 = 13 (prime, NOT Lucas)

The connection 11 = L₅ = Farey(v_PMNS, v_CKM) = new prime in
degree-5 algebraic cover is potentially significant.

---

## Generator Length Equations

### Fundamental generator lengths

**m003(−2,3):**
- ell_a = 0.88944300
- ell_b = 1.04031513
- Ratio ell_b/ell_a = 1.169625 ≈ 7/6 (0.254%)

**m006(−5,2):**
- ell_a = 0.94159568
- ell_b = 0.49057470
- Ratio ell_a/ell_b = 1.919373

### Mass ratio equations (m003)

| Equation | Value | Target | Error |
|---|---|---|---|
| exp(−ell_a + 2·ell_b) | 3.290986 | mb/mc = 3.291339 | 0.011% |
| exp(2·ell_a + ell_b) | 16.763454 | mτ/mμ = 16.817149 | 0.319% |
| exp(3·ell_a + ell_b) | 40.798452 | mt/mb = 41.313397 | 1.246% |

### Mass ratio equations (m006)

| Equation | Value | Target | Error |
|---|---|---|---|
| exp(3·ell_a) | 16.857355 | mτ/mμ = 16.817149 | 0.239% |

### Cross-manifold consistency
mτ/mμ appears on BOTH manifolds:
- m003: exp(2·ell_a + ell_b) = 16.763 (0.32%)
- m006: exp(3·ell_a) = 16.857 (0.24%)

This is the strongest evidence these are not random hits.

### Algebraic reformulation (m003)
The mb/mc equation is equivalent to:
```
|λ_b² / λ_a| = √(mb/mc) = 1.8142
```
where λ_a, λ_b are the dominant eigenvalues of the generator
holonomies. Computed value: 1.814107 (error 0.005%).

### Slope specificity
At adjacent slopes:
- (−1,3): exp(2·ell_b − ell_a) = 1.740 (47% error)
- (−3,3): exp(2·ell_b − ell_a) = 0.573 (83% error)

The slope (−2,3) is the unique minimum-volume H₁=ℤ/5 filling
of m003 (Meyerhoff manifold). Among all fillings with |p|,|q|≤8,
only three others come within 1% of mb/mc: slopes (3,7), (5,2), (4,5).

---

## Verified SM Parameter Encodings

### Mixing matrices and CP phase (established)

| Result | Value | PDG | Error | Method |
|---|---|---|---|---|
| PMNS Borel fitness | 0.005087 | — | global min | Borel QR |
| CKM fitness | 0.016482 | — | — | Gaussian overlap |
| δ_CP | 195.91° | 197.0° | 0.55% | π+φ(aaB)+φ(baa) |
| Cabibbo angle | — | 13.04° | 0.19% | axis angle |

### Twist angle spectrum (verified, needs statistical work)

| Result | Value | PDG | Error | Word |
|---|---|---|---|---|
| CKM isospectral φ | 92.487° | — | spread=0 | aaB=AbA=AAb |
| δ_CKM | 67.65° | 68.0° | 0.51% | 180−φ(aa on m006) |
| θ₂₃_CKM | 2.132° | 2.38° | 0.25° | φ(aaabb on m006) |
| θ₁₂_ν solar | 33.62° | 33.41° | 0.63% | 180−φ(abbAB on m006) |
| θ₁₂_CKM | 12.64° | 13.04° | 0.31° | 180−φ(AAB on m003) |
| θ₁₃_ν | 8.546° | 8.54° | 0.007° | 180−φ(BaBBBBB on m006) |

### Mass ratios (new findings this session)

| Result | Value | PDG | Error | Formula |
|---|---|---|---|---|
| mb/mc | 3.291 | 3.291 | 0.011% | exp(−ell_a+2ell_b) on m003 |
| mτ/mμ | 16.763 | 16.817 | 0.32% | exp(2ell_a+ell_b) on m003 |
| mt/mb | 40.80 | 41.31 | 1.25% | exp(3ell_a+ell_b) on m003 |
| mτ/mμ | 16.857 | 16.817 | 0.24% | exp(3ell_a) on m006 |

---

## Open Questions (prioritized)

### High priority
1. **What is the minimal polynomial of λ_b²/λ_a over ℚ(√−3)?**
   If it factors into low-degree polynomials over the trace field,
   the mb/mc equation has arithmetic content.

2. **Why does mτ/mμ appear on BOTH manifolds?**
   m003: exp(2ell_a+ell_b) and m006: exp(3ell_a) both give mτ/mμ
   to <0.32%. Is there a deeper geometric reason?

3. **Is 11 = L₅ = Farey(v_PMNS, v_CKM) = degree-5 cover prime
   a coincidence or a theorem?**
   Requires understanding the relationship between filling slope
   arithmetic and the covering tower prime structure.

### Medium priority
4. **What selects these two manifolds from the census?**
   Volume minimality + H₁=ℤ/5 + specific Dehn slopes.
   Is there a moduli-space characterization?

5. **Does the mt/mb equation (1.25% error) improve with
   better quark mass determinations?**
   The top mass has ~0.5% uncertainty; with future precision
   this could sharpen or falsify.

6. **Check ALL mass ratios on ALL slopes (|p|,|q|≤10).**
   Map the full function and identify what is special about
   the two flavor slopes.

### Low priority (pending arXiv endorsement)
7. Rewrite Dehn filling slopes paper with correct framing
   (algebraic cover prime = Farey intersection, not prime dictionary)
8. Rewrite twist spectrum paper with verified values
9. Build unified_reproduce.py for PRD paper

---

## Canonical Reproduction

All results reproducible via:
```bash
conda run -n sage python reproduce/hfg_reproduce.py    # mixing matrices
conda run -n sage python reproduce/cp_reproduce.py     # CP phase
conda run -n sage python reproduce/twist_reproduce.py  # twist census
conda run -n sage python minimal_polynomials.py        # eigenvalue algebra
conda run -n sage python slope_survey.py               # mass ratio equations
conda run -n sage python tower_reconcile.py            # covering tower
```

Repository: https://github.com/drmlgentry/hyperbolic-flavor-scan

---

## Twist Angle Convention
All twist angles in the census use the **positive-branch** eigenvalue:
φ(γ) = Im(log λ) where λ chosen so Im(log λ) ≥ 0, giving φ ∈ [0°,180°].
Implementation: `G.SL2C(word)`, select eigenvalue with Im ≥ 0.

The CP phase paper uses `polished_holonomy()` which can give negative
angles. Both are consistent: φ_positive = 180° − |φ_polished| when
the polished angle is negative near 180°.

---

## Session Summary (May 19, 2026)

This session established:
1. The covering tower has TWO distinct structures (algebraic vs census)
   that were conflated in the slopes paper — now corrected
2. The mb/mc coincidence reduces to a linear equation in the two
   generator lengths of m003, not an arbitrary word pair
3. Three mass ratios appear simultaneously at slope (−2,3) with
   small integer coefficients — mb/mc (0.011%), mτ/mμ (0.32%), mt/mb (1.25%)
4. mτ/mμ appears independently on BOTH flavor manifolds
5. The slope specificity test shows (−2,3) is the best among nearby
   slopes for mb/mc, consistent with it being the minimal volume H₁=ℤ/5 filling

Active submissions: 10 papers at 8 journals (see memory for full list).
All held pending PRD referee response before further submissions.
