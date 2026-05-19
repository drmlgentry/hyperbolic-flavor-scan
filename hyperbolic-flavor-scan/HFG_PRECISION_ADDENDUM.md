# HFG Precision Analysis — Addendum to May 2026 Elucidation

## High-Precision Test Results (May 19, 2026)

### Calibration

The SnapPy `polished_holonomy()` computation achieves approximately
**16 bits** of reliable precision for algebraic numbers, even when
run at 1000-bit internal arithmetic. This was established by comparing
the computed Re(tr(ρ(aa) on m006)) to the known exact value 3−√17:

```
Computed:  -1.12311843562865...
Exact:     -1.12310562561766...
Error:      1.28e-5  (~16 bits of agreement)
```

This 16-bit ceiling applies to all holonomy-derived quantities.

### Verdict on Mass Ratio Claims

| Claim | Error | Bits | Verdict |
|---|---|---|---|
| \|λ_b²/λ_a\|² = mb/mc | 0.011% | 13.2 | **INCONCLUSIVE** — within calibration ceiling |
| exp(3ℓ_a) = mτ/mμ (m006) | 0.239% | 8.7 | **LIKELY APPROXIMATE** — 7 bits below calibration |
| exp(2ℓ_a+ℓ_b) = mτ/mμ (m003) | 0.319% | 8.3 | **LIKELY APPROXIMATE** — 7 bits below calibration |
| exp(3ℓ_a+ℓ_b) = mt/mb (m003) | 1.246% | 6.6 | **APPROXIMATE** — clearly coincidental |

### Interpretation

**mb/mc (inconclusive):** The 13.2-bit agreement is within the 16-bit
calibration ceiling. We cannot distinguish an exact relation from a
numerical coincidence at this precision level. The result is consistent
with either interpretation.

**mτ/mμ (likely approximate):** Both manifolds give ~8-9 bits of
agreement, significantly below the 16-bit calibration. A genuinely
exact algebraic relation would agree to at least calibration precision.
The cross-manifold "consistency" appears to be coincidental matching
at the ~0.3% level.

**mt/mb (approximate):** 1.25% error rules out any exact relation.

### What This Changes

The mass ratio equations from generator lengths are:
- **Not demonstrated to be structurally exact**
- May be numerical coincidences at the 0.01%–1.25% level
- Should NOT be claimed as predictions in papers

The strong published results remain unaffected:
- CKM isospectrality (theorem, exact)
- CP phase 195.91° (verified to holonomy precision)
- PMNS Borel fitness 0.005087 (global minimum, p<10⁻⁴)
- Covering tower prime 11 = Farey intersection = L₅ (algebraically exact)

### Path to Resolution

**Option A — Higher precision holonomy:**
SnapPy's `snap()` method uses exact arithmetic in the trace field and
may achieve higher precision than `polished_holonomy()`. Computing
tr(ρ(γ)) exactly via the trace field arithmetic would give a
definitive answer.

**Option B — Theoretical argument:**
Derive from hyperbolic Dehn surgery theory whether the generator
lengths of m003(−2,3) satisfy the equation
2ℓ_b − ℓ_a = ln(mb/mc) as an exact consequence of the slope
arithmetic, or only approximately.

**Option C — Comparative study:**
Check whether nearby fillings (e.g., m003(−3,5), m003(−3,7)) also
have generator lengths encoding mass ratios at similar precision.
If yes: generic, not special. If no: specific to (−2,3).

---

## Corrected Working Summary

### Confirmed structural results:
1. **CKM axis isospectrality** — theorem from H₁ class collision [AbA]=[AAb]=2
2. **CP phase δ=195.91°** — zero free parameters, 0.55% from PDG, structurally motivated
3. **PMNS global minimum fitness** — 0.005087, p<10⁻⁴ vs Haar random
4. **Covering tower** — unique degree-5 cover has H₁=ℤ/11+ℤ/11, prime 11=L₅=Farey intersection

### Numerologically suggestive but unconfirmed:
5. **mb/mc from generator lengths** — 0.011%, inconclusive at current precision
6. **mτ/mμ cross-manifold** — 0.24-0.32%, likely approximate

### Disconfirmed:
7. **mt/mb** — 1.25%, coincidental
8. **MZ/MW statistical significance** — does not survive proper null test
9. **"Prime dictionary" from slopes** — conflates algebraic covers with census coincidences
