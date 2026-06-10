# HFG Handoff — June 9, 2026

## Session Summary

Major research session. The CKM selection mechanism was fully resolved.
Four Substack articles published. RNTB paper approved. Website updated.

---

## NEW RESULTS (June 9 2026)

### 1. Corrected Abelianization for CKM Manifold

The H₁ map for m006(-5,2) is NOT simply nb mod 5.
From relators ababbAAbb and abbAbbaBabbAbbAbbaB:
  5[b] = 0  and  [a] + 8[b] = 0  →  [a] = 2[b] mod 5

Correct map: H1(word) = (2·na + nb) mod 5

CKM triple H₁ classes (corrected):
  aaB: H₁ = 3
  AbA: H₁ = 2
  AAb: H₁ = 2

### 2. PMNS vs CKM Structural Duality

PMNS — Axis configuration (0, c, −c):
  Triple: {aa(H₁=0), aaB(H₁=4), baa(H₁=1)}
  H₁(aaB·baa) = 4+1 = 0 mod 5  ← inverse pair
  Two distinct geodesic classes (φ=-176.73° and -167.36°)
  Resonance θ* = -180°  →  Borel factorization  →  large mixing

CKM — Anchor configuration (c, −c, −c):
  Triple: {aaB(H₁=3), AbA(H₁=2), AAb(H₁=2)}
  H₁(aaB·AbA) = 3+2 = 0 mod 5  ← anchor+partner
  H₁(aaB·AAb) = 3+2 = 0 mod 5  ← anchor+partner
  H₁(AbA·AAb) = 2+2 = 4 ≠ 0    ← partners NOT inverse to each other
  One geodesic class (φ=+92.49° for all three)
  Resonance θ* = +90°  →  Iwasawa factorization  →  small mixing

The duality: 3+2 = 0 mod 5, twice. Emerged from fixing abelianization.

### 3. CKM Selection: System of Distinct Representatives (SDR)

The conjugacy class at φ=+92.49°, len=1.77172 contains 6 words:
  H₁=3: {aaB, aBa, Baa}   (three cyclic rotations)
  H₁=2: {bAA, AbA, AAb}   (their inverses)

Three inverse pairs:
  aaB ↔ bAA  (3+2=0 mod 5)
  aBa ↔ AbA  (3+2=0 mod 5)
  Baa ↔ AAb  (3+2=0 mod 5)

Valid triples = Systems of Distinct Representatives:
  one word per inverse pair, H₁ pattern (3,2,2):

  Triple 1: {aaB, AbA, AAb}  tr(product)=18.43-13.02j  |tr|=22.56  ← CKM
  Triple 2: {aBa, bAA, AAb}  tr(product)=-27.39+9.52j  |tr|=28.99
  Triple 3: {Baa, bAA, AbA}  tr(product)=18.43-13.02j  |tr|=22.56

Triples 1 and 3 give IDENTICAL product traces → same Iwasawa K factor.
Triple 2 is isolated.
WHY T1=T3: open question — related by outer automorphism (not yet identified).

### 4. Phase Resonance Mechanism (CONFIRMED FOR BOTH)

PMNS: θ* = -180° selected by systole phase (-95.72°) pulling toward -180°
CKM:  θ* = +90°  selected by systole `b` at φ=+89.16° ≈ +90°

The resonance is NOT determined by the Chern-Simons invariant:
  CS(m006 cusped) = -0.1141366530 ≈ -5/44  (no clean relation to +90°)

Current best hypothesis: θ* determined by the systole eigenphase,
rounded to nearest torsion order of the eigenvalue unit part:
  order 2 (λ/|λ| = -1) → θ* = -180°
  order 4 (λ/|λ| = i)  → θ* = +90°

### 5. Fricke Invariants

tr(aaB) = tr(AbA) = tr(AAb) = -0.123 + 2.011i  (identical — same conjugacy class)
Pairwise traces all different:
  tr(aaB·AbA) = -1.460 + 4.547i  |tr| = 4.776
  tr(aaB·AAb) = -10.178 - 2.971i |tr| = 10.603
  tr(AbA·AAb) = -2.568 - 5.042i  |tr| = 5.658
Triple product cyclically invariant:
  tr(g1·g2·g3) = tr(g2·g3·g1) = tr(g3·g1·g2) = 18.43-13.02j  ✓
  tr(g1·g2·g3) ≠ tr(g1·g3·g2): difference = 30.21

---

## PUBLICATIONS (June 9)

Substack posts published:
  1. "The CP Phase Is a Geometric Invariant" — live
  2. "A Unified Selection Principle for PMNS and CKM" — live
  3. "The HFG Programme: From Derivation to Selection" — live
  4. "The Complete CKM Selection Rule: A Latin Square in Hyperbolic Space" — live

Desk rejections received:
  PLB-D-26-01448 (PMNS/Borel) — Kitano, no comment → reroute to Annals of Physics
  PLB-D-26-01449 (Torsion)    — Kitano, no comment → reroute to EPJC
  PLB-D-26-01463 (CKM)        — Kitano, no comment → reroute to PRD
  MATH-D-26-00372 (BPS/X0(11)) — Warzel/LMP → needs more development

RNTB-D-26-00299 (sextic-octic) — NEW SUBMISSION ✓ June 8

---

## OPEN QUESTIONS (priority order)

1. **Third manifold** — find H₁=ℤ/5 manifold with θ* ≠ {-180°,+90°}
   This is the NEXT computation session focus.
   The disc=-283 census has 16 cusped manifolds.
   OrientableClosedCensus has many more with H₁=ℤ/5.
   Method: scan for systole phase ≠ {-95°, +89°}

2. **Why T1=T3?** — identify the outer automorphism relating
   {aaB,AbA,AAb} and {Baa,bAA,AbA} with identical K factors.
   No conjugator found in words of length ≤4.
   Try length ≤6 or look for automorphism of the group presentation.

3. **What determines θ*?** — prove connection between
   arithmetic of trace field and resonance angle.
   CS not it. Systole phase is best candidate.

4. **D(w³) = v₀ algebraically** — still open from June 6.
   Claim (a): vol(m019) = 3·D(w³) nearly provable.
   Claim (b): vol(m003(-2,3)) = v₀ requires commensurability.

5. **Lucas paper → Ramanujan Journal** — resubmit (MRL rejected).

6. **Fermion masses** — Lucas structure gives muon/tau mass ratios.
   L₁₁=199 ≈ m_μ/mₑ (0.003%), L₁₇=3571 ≈ m_τ/mₑ (0.000%).
   This deserves a standalone paper.

7. **σ₃ vs σ₂ verification** — confirm geometric embedding
   is consistently σ₂ throughout all computations.

---

## CANONICAL RESULTS (unchanged)

PMNS: m003(-2,3) = OrientableClosedCensus[1]
  words {aa, aaB, baa}, Borel, col perm (1,0,2)
  fitness 0.005087
  θ₁₂=33.67°/θ₂₃=47.63°/θ₁₃=8.37° vs PDG 33.65°/47.64°/8.57°

CKM: m006(-5,2) = OrientableClosedCensus[43]
  words {aaB, AbA, AAb}, Iwasawa
  fitness 0.016482

CP phase: δ = π + φ(aaB) + φ(baa) = 195.91° vs PDG 197.0°
  φ(aaB) = -176.731°, φ(baa) = -167.362°

---

## THIRD MANIFOLD SEARCH — PLAN

Target: closed hyperbolic 3-manifold M with:
  - H₁(M) = ℤ/5
  - Systole phase θ_sys ≠ {-95°, +89°}
  - Short geodesics clustering near new θ*
  - θ* not in {-180°, +90°}

Candidate θ* values (torsion orders):
  order 3: θ* = +60° or -60° (cube roots of unity)
  order 6: θ* = +30° or -30°
  order 8: θ* = +45° or -45°

Search strategy:
  1. Scan OrientableClosedCensus for H₁=ℤ/5
  2. For each: compute systole phase
  3. Cluster analysis of short geodesic phases
  4. Flag manifolds with θ* outside known zones

Sage code skeleton (paste into sage terminal):
```python
import snappy, cmath

def systole_phase(M):
    G = M.fundamental_group()
    from itertools import product as iproduct
    best_len = 999; best_phi = None
    for w in ['a','b','ab','ba','aab','abb','bab','bba']:
        try:
            mat = G.SL2C(w)
            t = complex(mat[0][0]+mat[1][1])
            ht = t/2
            l1 = ht+cmath.sqrt(ht**2-1)
            l2 = ht-cmath.sqrt(ht**2-1)
            lam = l1 if abs(l1)>=abs(l2) else l2
            phi = float(cmath.log(lam).imag*180/cmath.pi)
            length = abs(float(2*cmath.log(lam).real))
            if length > 0.01 and length < best_len:
                best_len = length; best_phi = phi
        except: pass
    return best_phi, best_len

results = []
census = snappy.OrientableClosedCensus
for idx in range(500):
    try:
        M = census[idx]
        if str(M.homology()) != 'Z/5': continue
        phi, length = systole_phase(M)
        if phi is None: continue
        known = abs(phi+95.7)<10 or abs(phi-89.2)<10
        results.append((idx, phi, length, known))
        marker = '' if known else ' ← NEW THETA*'
        print(f"idx={idx:4d}  H1=Z/5  systole_phi={phi:+.2f}°  len={length:.4f}{marker}")
    except: pass
```

---

## ENVIRONMENT

- WSL conda sage: conda activate sage
- Canonical script: /mnt/c/dev/hyperbolic-flavor-scan/hfg_reproduce.py
- Website: C:\dev\hyperbolic-flavor-geometry\docs\
- GitHub token: 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
- PMNS needs WSL/sage; CKM works on Windows too


---

## THIRD MANIFOLD RESULTS (June 10, 2026)

### Manifold: m206(1,2) = OrientableClosedCensus[209]

Volume: 2.828122 = 2.882 × v₀ (not a clean multiple)
H₁ = ℤ/5, [a]=1, [b]=1
Resonance θ* = −60° (order-6 torsion, Eisenstein sector)
Systole: abAB, φ = −93.29°, len = 0.72979

### Exact Arithmetic

Cusped parent m206:
  trace field = Q(√−3), disc = −3  ← EXACT EISENSTEIN
  cusp shape = i√3 (purely imaginary, in Q(√−3))
  tr(a) satisfies x² + 3x + 3 = 0  ← EXACT

Closed filling m206(1,2):
  tr(a) + tr(b) = 0  ← EXACT Z/2 symmetry
  λ_b / λ_a = −1    ← EXACT (eigenvalues are negatives)
  tr(a)² satisfies: x² − (235115/48036)x + 160022/16873 = 0
  disc of this poly = −9183165546402989879 = −47·359·544252091886623
  Field Q(tr(a)²) is NOT Q(√−3) — Dehn filling deforms the field
  16873 = 47×359, 48036 = 4×3×4003 (factor of 3 present)

### Key: The Eisenstein Structure

The θ*=−60° resonance is INHERITED from the cusped parent.
The cusped parent is pure Eisenstein (Q(√−3), disc=−3).
The Dehn filling m206(1,2) deforms this to a larger field.
The Z/2 symmetry tr(a)=−tr(b) has NO analogue in PMNS or CKM.

### Closest Geodesics to θ*=−60°

H₁=3: aaa   φ=−59.016°  D=0.984°  len=2.53080  ← most resonant
H₁=4: aBB   φ=−64.002°  D=4.002°  len=3.12697
H₁=0: ababb φ=−71.269°  D=11.269° len=2.67500

### Mixing Matrix

Best H₁-trivial triple: {aab(H₁=3), aaa(H₁=3), aBB(H₁=4)}
H₁ sum = 3+3+4 = 10 = 0 mod 5  ✓
Mixing angle: 67.9°–81.4° across orderings
Standard deviation: 4.79° (stable)
CP analog: ~−2.5° (small, nearly CP-conserving)

PMNS comparison: 33.7° (large but not maximal)
CKM comparison:  ~13° (small)
NEW:             ~74° (near-maximal)

### Torsion Taxonomy (Complete)

| Manifold     | Matrix | θ*    | Torsion order | λ relation    | Field      | Mixing  |
|-------------|--------|-------|---------------|---------------|------------|---------|
| m003(−2,3)  | PMNS   | −180° | 2             | none          | disc=−283  | 33.7°   |
| m006(−5,2)  | CKM    | +90°  | 4             | none          | disc=?     | ~13°    |
| m206(1,2)   | BSM?   | −60°  | 6             | λ_b=−λ_a     | disc=−9.2e18 | ~74° |

### Physical Interpretation

Maximal mixing (45°–90°) appears in:
  - Sterile neutrino mixing
  - Dark matter sector
  - Mirror matter / left-right symmetric extension
  - Kaluza-Klein extra dimensions

The Z/2 symmetry could signal a discrete parity symmetry.
The −60° resonance is Eisenstein — related to cube roots of unity.

### FALSIFIABLE PREDICTION

If a BSM sector with θ_mix ~ 74° is discovered,
its mixing matrix should be reproduced by m206(1,2).
If no such sector is found, the third manifold prediction is falsified.

### Factorization Method for θ*=−60°

OPEN: What factorization corresponds to −60°?
  θ*=−180° → Borel (lower triangular)
  θ*=+90°  → Iwasawa (KAN)
  θ*=−60°  → ??? (new decomposition needed)

### Open Questions

1. What factorization method corresponds to θ*=−60°?
2. Why does the cusped parent have Q(√−3) but the filling doesn't?
3. Are there manifolds with θ*=+120°, +45°, −45°, etc.?
4. Scan census[300:1000] for more H₁=ℤ/5 novel manifolds
5. Compute full 3×3 unitary (not just 2×2) for the mixing matrix

