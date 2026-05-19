# HFG Research Program — Master Handoff Document
**Author:** Marvin L. Gentry  
**Last updated:** May 19, 2026  
**Repository:** https://github.com/drmlgentry/hyperbolic-flavor-scan

---

## What This Program Is

The Hyperbolic Flavor Geometry (HFG) program proposes that Standard
Model flavor parameters arise from the arithmetic geometry of two
specific compact hyperbolic 3-manifolds. This document is the
authoritative record of what has been established, what is open,
and how to reproduce every result.

---

## The Two Manifolds — Complete Portrait

### M_PMNS = m003(−2,3)
```
SnapPy:       OrientableClosedCensus[1]   ← ALWAYS USE INDEX
Volume:       0.98136883
H₁:           ℤ/5
Tetrahedra:   2
Orientable:   True
Chern-Simons: 1/4
Relators:     ababAbbAb, abAbaabAbaBAB
Trace field:  ℚ(√−3)  [imaginary quadratic]
  Ring Z[ω], ω = e^{2πi/3}, norm a²−ab+b²
  Discriminant −3 < 0 → loxodromic holonomy → large CP violation
Dehn slope:   (−2,3) on cusped manifold m003
```

### M_CKM = m006(−5,2)
```
SnapPy:       OrientableClosedCensus[43]  ← ALWAYS USE INDEX
Volume:       2.02885309
H₁:           ℤ/5
Tetrahedra:   3
Orientable:   True
Relators:     ababbAAbb, abbAbbaBabbAbbAbbaB
Trace field:  ℚ(√17)  [real quadratic]
  Ring Z[√17], norm a²−17b²
  Discriminant +17 > 0 → real/hyperbolic holonomy → CP suppression
  Verified: tr(ρ(aa)) = 3−√17 (MSbar)
Dehn slope:   (−5,2) on cusped manifold m006
```

### Slope Arithmetic
```
Filling slopes: v_PMNS = (−2,3),  v_CKM = (−5,2)
Farey intersection: det[(−2,−5),(3,2)] = −4+15 = 11 = L₅ (Lucas prime)
PMNS slope norm: ‖(−2,3)‖² = 4+9 = 13  (prime, NOT Lucas)
CKM slope norm:  ‖(−5,2)‖² = 25+4 = 29 = L₇ (Lucas prime)
```

---

## Generator Holonomy Data

### m003(−2,3) — verified at 150-bit via polished_holonomy()

| word | ℓ | φ (°) | \|λ\| |
|---|---|---|---|
| a   | 0.88944 | 84.278  | 1.56006 |
| b   | 1.04032 | 28.143  | 1.68229 |
| aa  | 1.77889 | 168.556 | 2.43377 |
| bb  | 2.08063 | 56.286  | 2.83011 |
| aaB | 1.73425 | −176.731| 2.38006 |
| baa | 1.99665 | −167.362| 2.71373 |
| AbA | 1.73425 | −176.731| 2.38006 |  ← isospectral to aaB
| AAb | 1.73425 | −176.731| 2.38006 |  ← isospectral to aaB
| bAbA| 1.77889 | 168.556 | 2.43377 |  ← isospectral to aa

### m006(−5,2) — verified at 150-bit

| word  | ℓ | φ (°) | \|λ\| |
|---|---|---|---|
| a    | 0.94160 | 56.173  | 1.60127 |
| b    | 0.49057 | −90.844 | 1.27799 |
| aa   | 1.88319 | 112.346 | 2.56407 |
| aaB  | 1.77172 | −87.513 | 2.42506 |
| AbA  | 1.77172 | −87.513 | 2.42506 |  ← isospectral to aaB
| AAb  | 1.77172 | −87.513 | 2.42506 |  ← isospectral to aaB
| aaabb| 3.46544 | 2.132   | 5.65601 |
| abbAB| 1.72900 | 33.620  | 2.37382 |

**CKM isospectrality:** aaB = AbA = AAb have identical complex lengths.
Spread = 0.000000° exactly. This is a theorem from H₁ class collision
[AbA] = [AAb] = 2 in H₁ = ℤ/5, forcing J_CKM = 0.

---

## Twist Angle Convention

**Positive-branch convention** (used in twist census):
```
φ(γ) = Im(log λ)  where λ chosen so Im(log λ) ≥ 0
φ ∈ [0°, 180°]
Implementation: G.SL2C(word), select eigenvalue with Im ≥ 0
```

**Polished holonomy convention** (used in CP phase):
```
φ(γ) = Im(log λ)  where λ is dominant eigenvalue |λ| > 1
φ can be negative, ∈ (−180°, 180°]
Implementation: M.polished_holonomy() at 150-bit
```

Both conventions are consistent: φ_positive = 180° − |φ_polished|
when the polished angle is negative near −180°.

---

## Verified Results — Strength Hierarchy

### TIER 1: Exact theorems (algebraically verified)

**T1.1 CKM isospectrality**
All three CKM axis words have identical complex lengths on m006:
- ℓ(aaB) = ℓ(AbA) = ℓ(AAb) = 1.77172, φ = −87.513° (spread = 0)
- Proof: H₁ class collision [AbA] = [AAb] = 2 mod 5 forces conjugacy
- J_CKM = 0 follows geometrically

**T1.2 Covering tower prime**
M_PMNS has a unique degree-5 algebraic cover with H₁ = ℤ/11 + ℤ/11.
New prime introduced: 11 = L₅.
Verified by M_PMNS.covers(5) returning exactly one cover.
Farey intersection of the two filling slopes = 11 = L₅ (same prime).

**T1.3 Power law**
ℓ(γⁿ) = n·ℓ(γ) and φ(γⁿ) = n·φ(γ) mod 360° exactly.
All integer ratios in the length spectrum follow from this.

**T1.4 Trace field verification**
tr(ρ(aa) on m006) = 3−√17 verified to 16-bit holonomy precision.
Satisfies x² − 6x − 8 = 0 (roots 3 ± √17).

---

### TIER 2: Computationally verified, structurally motivated

**T2.1 CP phase**
δ_CP = π + φ(aaB) + φ(baa) mod 360° = 195.91°
PDG 2024: 197.0°, error 0.55%, zero free parameters.
Sign fixed geometrically: det[n̂(aa), n̂(aaB), n̂(baa)] = +0.079 > 0.
Formula: `(180 + phi_aaB + phi_baa) % 360` using polished holonomy.

**T2.2 PMNS Borel fitness**
Borel (lower-triangular N-factor) construction on m003:
- Word triple {aa, aaB, baa}: fitness 0.005087
- Global minimum: verified over 50,000 Haar-random unitaries (p < 10⁻⁴)
- Symmetric QR floor: 0.300 (Borel is required, not optional)

**T2.3 CKM Gaussian fitness**
Gaussian overlap (K-factor) on m006:
- Word triple {aaB, AbA, AAb}: fitness 0.016482, σ = 0.49
- Cabibbo angle error: 0.19%

**T2.4 Twist angle census coincidences**
Verified values (positive-branch convention on m006 unless noted):
| Observable | Word | Value | PDG | Error |
|---|---|---|---|---|
| δ_CKM | aa (m006) | 180°−φ = 67.65° | 68.0° | 0.51% |
| θ₂₃_CKM | aaabb (m006) | φ = 2.132° | 2.38° | 0.25° |
| θ₁₂_ν solar | abbAB (m006) | 180°−φ = 33.62° | 33.41° | 0.63% |
| θ₁₂_CKM | AAB (m003) | 180°−φ = 12.64° | 13.04° | 0.31° |
| θ₁₃_ν | BaBBBBB (m006) | 180°−φ = 8.546° | 8.54° | 0.007° |

Statistical note: census enrichment factor ~1.27x over random.
Individual matches are suggestive but not statistically significant
without geometric selection principle (folding is ad hoc without
orientation justification as in T2.1).

---

### TIER 3: Numerically striking, not yet structurally explained

**T3.1 mb/mc from generator lengths (m003)**
|λ_b²/λ_a|² = 3.290986 vs PDG mb/mc = 3.291339
Error: 0.011%, within 0.006σ of experimental central value.

Key facts:
- Reduces to: 2ℓ_b − ℓ_a = ln(mb/mc) at 0.011%
- Algebraically: |λ_b²/λ_a|² = exp(2ℓ_b − ℓ_a)
- Slope-specific: adjacent slopes give 47-83% errors
- High-precision test (1000-bit): 13.2 bits agreement
- Exactness: NOT RULED OUT (within experimental mb/mc uncertainty)
- Status: striking coincidence demanding geometric explanation

**T3.2 Generator length ratio**
ℓ_b/ℓ_a = 1.04032/0.88944 = 1.16963 on m003
≈ 7/6 to 0.25% — no geometric explanation yet

---

### TIER 4: Disconfirmed or ruled out

**T4.1 mτ/mμ cross-manifold** — 0.24-0.32% error, 8-9 bits precision
Significantly below calibration precision (16 bits). Likely approximate.

**T4.2 mt/mb from generator lengths** — 1.25% error. Coincidental.

**T4.3 MZ/MW statistical significance** — does not survive null test
on distinct geodesic lengths. p = 0.32 (not significant).

**T4.4 "Prime dictionary" from slope arithmetic** — conflates
algebraic covering spaces with census volume coincidences.
Census manifolds at n·vol are commensurable neighbors, NOT covers.
Algebraic covers through degree 6: only one (degree 5, prime 11).

---

## Covering Tower — Corrected Picture

### Method A: Algebraic (M.covers())
Finds ALL finite-sheeted covering spaces algebraically.
Through degree 6, M_PMNS has ONE cover:
- Degree 5: H₁ = ℤ/11 + ℤ/11, new prime = 11 = L₅

M_PMNS is arithmetically rigid.

### Method B: Census volume matching
Finds census manifolds with volume = n·vol(M_PMNS).
These are commensurable neighbors, NOT algebraic covers.
- Degree 2: H₁ primes {5,11} and {7} — all Lucas primes
- Degree 3: H₁ primes {2,3,5,7} — all Lucas primes

### The genuine connection
- Algebraic cover prime 11 = L₅
- Farey intersection of slopes = det[(−2,−5),(3,2)] = 11 = L₅
- CKM slope norm = 29 = L₇
This is the correct framing for the covering tower paper.

---

## Canonical Reproduction Commands

```bash
# Environment
conda activate sage
cd C:\dev\hyperbolic-flavor-scan

# PMNS Borel fitness (0.005087) and CKM Gaussian fitness (0.016482)
conda run -n sage python reproduce/hfg_reproduce.py

# CP phase (195.91°, zero free parameters)
conda run -n sage python reproduce/cp_reproduce.py

# Twist angle census (9 claims, PASS/FAIL table, ~3 min)
conda run -n sage python reproduce/twist_reproduce.py

# Manifold elucidation (geometry, spectrum, ratios)
conda run -n sage python elucidate_simple.py

# Tower reconciliation (algebraic vs census covers)
conda run -n sage python tower_reconcile.py

# Slope survey (mb/mc across Dehn fillings)
conda run -n sage python slope_survey.py

# High-precision test (mb/mc exactness, 1000-bit)
conda run -n sage sage high_precision_test.sage

# Algebraic number recognition
conda run -n sage sage algdep_focused.sage
```

Expected runtimes: hfg_reproduce ~2 min, twist_reproduce ~3 min,
slope_survey ~5 min, high_precision_test ~3 min.

---

## Active Submission Portfolio (May 19, 2026)

| ID | Journal | Paper | Status |
|---|---|---|---|
| DS14327 | PRD | Unified HFG | With editors |
| PLB-D-26-01341 | PLB | CP phase | With editor |
| PLB-D-26-01006 | PLB | Charge conjugation | Under review |
| AHPO-D-26-00255 | AHP | PMNS mixing | New submission |
| AHPO-D-26-00231 | AHP | CP holonomy | New submission |
| JGP13076 | JGP | CKM mixing | Submitted |
| JPhysA-124843 | JPA | Discrete mixing operators | Submitted |
| TRGR → rejected | — | Discrete mixing (old) | Rejected |
| JMP26-AR-01272 | JMP | Qubit gates | 17+ days |
| JKTR | JKTR | Alexander polynomial | Submitted |

**Strategy:** Hold all new submissions pending PRD response.
Use waiting time to complete manifold description and fix papers
with verified values.

---

## Papers Needing Updates Before Next Submission

| Paper | Issue | Action |
|---|---|---|
| Twist spectrum (RINP-D-26-00330) | Old CP formula (203.5°), wrong MZ/MW words, wrong θ₁₃ source | Rewrite with verified values |
| Dehn filling slopes (SSRN 6761981) | Conflates algebraic covers with census coincidences | Reframe: algebraic cover prime = Farey intersection = L₅ |
| Lucas structure (SSRN 6754501, removed) | Tower primes {2,3,7,11,29} → {2,3,11} | Correct and repost |

---

## Open Questions (prioritized)

1. **mb/mc exactness** — Requires either SnapPy `snap()` exact
   arithmetic or a theoretical derivation from slope (−2,3) geometry.
   Currently: 0.011% match, within experimental uncertainty, not ruled out.

2. **Why 11 = Farey = cover prime?** — Is there a theorem connecting
   the Farey intersection of filling slopes to the first new prime in
   the algebraic covering tower? Pure mathematics question.

3. **Why these two manifolds?** — Minimal volume + H₁=ℤ/5 + specific
   slopes. Is there a moduli-space characterization that selects them?

4. **Folding justification** — The twist angle census uses 180°−φ
   without systematic geometric justification (except for the CP phase,
   where the determinant argument fixes the sign). A geometric selection
   principle for the census folding would strengthen Paper B.

5. **Algebraic covers at degree 7-19** — tower_extended.py found
   primes {2,3,11} through degree 19. Do primes 7 and 29 ever appear?
   If not, the Lucas-purity claim is strong through all computed degrees.

---

## Key File Locations

```
C:\dev\hyperbolic-flavor-scan\
  reproduce/
    hfg_reproduce.py     ← CKM + PMNS fitness (canonical)
    cp_reproduce.py      ← CP phase 195.91°
    twist_reproduce.py   ← twist census 9 claims
  
  HFG_ELUCIDATION_MAY2026.md  ← geometric findings
  HFG_PRECISION_ADDENDUM.md   ← precision analysis
  HANDOFF.md                  ← this document
  REPRODUCE.md                ← referee instructions
  
  slope_survey.py         ← mass ratios across Dehn fillings
  tower_reconcile.py      ← algebraic vs census covers
  elucidate_simple.py     ← manifold portrait
  high_precision_test.sage← mb/mc exactness test
  algdep_focused.sage     ← algebraic number recognition

C:\dev\framework\papers\
  hfg-unified\            ← PRD submission (DS14327)
  hyperbolic-flavor-ckm-v2\  ← JGP submission (JGP13076)
  hyperbolic-flavor-pmns\    ← AHP submission (AHPO-D-26-00255)
  hyperbolic-flavor-cp-v2\   ← PLB submission (PLB-D-26-01341)
  discrete-mixing-jpa\       ← JPA submission (JPhysA-124843)
```

---

## Session Log Summary

| Date | Key findings |
|---|---|
| Mar 2026 | Initial four RINP submissions; m003/m006 identified |
| Apr-May 2026 | All four HFG pillars established: CKM, PMNS, CP, mass norms |
| May 16 | PRD unified paper submitted (DS14327) |
| May 17 | PLB CP phase submitted; AHP PMNS submitted; JGP CKM submitted |
| May 18 | SSRN updates; TRGR rejection; JPA discrete mixing submitted |
| May 19 | Deep geometric elucidation: generator length equations, tower reconciliation, precision analysis. mb/mc: 0.011%, slope-specific, exactness inconclusive. |
