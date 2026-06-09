# hyperbolic-flavor-scan

**Scan code for the HFG programme — SnapPy/Sage census scans, holonomy phase analysis, geodesic enumeration, Lucas trace, Galois-Weyl correspondence (2026).**

Part of [Hyperbolic Flavor Geometry](https://hyperbolicflavorgeometry.org) by Marvin L. Gentry.

📦 **[LatticeFit v0.2.0 on PyPI](https://pypi.org/project/latticefit/0.2.0/)**

---

## Quick Start

```bash
pip install snappy latticefit numpy scipy
python hfg_reproduce.py
```

**Expected output:**
```
PMNS: manifold m003(-2,3), fitness 0.005087  ✓
CKM:  manifold m006(-5,2), fitness 0.016482  ✓
CP phase: 195.91° (PDG 197.0°, error 1.09°)  ✓
```

**Platform**: PMNS requires WSL + `conda activate sage` (SnapPy 3.3.2). CKM reproduces on Windows and Linux. Always use `OrientableClosedCensus` by index, never by name string.

---

## Key Scripts

| Script | Purpose |
|---|---|
| `hfg_reproduce.py` | Canonical reproduction: PMNS + CKM fitness, CP phase |
| `census_scan.py` | Full SnapPy census scan for fitness landscape |
| `geodesic_phase_scan.py` | Eigenphase enumeration by H₁ class |
| `lucas_trace_check.py` | Verify Lucas prime structure in length spectrum |
| `galois_weyl_check.py` | Verify Gal(m003)=ℤ/2, Gal(m006)=S₃, Gal(m019)=S₄ |

---

## Canonical Results

### PMNS Manifold: m003(−2,3) = OrientableClosedCensus[1]

```python
M = snappy.ManifoldHP("m003(-2,3)")
# isosig: cPcbbbdxm_BaBb(-2,3)
# words: {aa, aaB, baa}
# factorization: Borel, column permutation (1,0,2)
# fitness: 0.005087
# θ₁₂=33.67° / θ₂₃=47.63° / θ₁₃=8.37°  vs PDG 33.65° / 47.64° / 8.57°
```

### CKM Manifold: m006(−5,2) = OrientableClosedCensus[43]

```python
M = snappy.ManifoldHP("m006(-5,2)")
# isosig: dLQacccjnjs_aBbB(-5,2)
# words: {aaB, AbA, AAb}
# factorization: Iwasawa
# fitness: 0.016482
```

### CP Phase (zero free parameters)

```
δ = π + φ_dom(aaB) + φ_dom(baa) = 195.91°
φ_dom(aaB) = −176.731°  (dominant eigenvalue branch, |λ| > 1)
φ_dom(baa) = −167.362°
PDG: 197.0°  error: 1.09%
```

---

## Selection Principle (June 2026)

The CP phase is a **manifold invariant**. The ℤ/5 phase resonance mechanism:

```
Phase resonance θ* = −180°  (real negative trace)
H₁ classes and D = |φ + 180°|:
  class 1: aab,   φ = −167.36°,  D = 12.64°  ← PMNS
  class 2: aaBABB, φ = −159.19°, D = 20.81°
  class 3: aaBAB,  φ = −151.86°, D = 28.14°
  class 4: aaB,   φ = −176.73°,  D =  3.27°  ← PMNS

Inverse pair (1,4): D-sum = 15.91°
Inverse pair (2,3): D-sum = 48.95°
Ratio: 3.1× — selected WITHOUT reference to PDG
```

---

## H₁ Structure

```
H₁(m003(−2,3)) = ℤ/5
[a] = 0  (a is in the commutator subgroup)
[b] generates ℤ/5
aa  → class 0 (axis word)
aaB → class 4 (CP phase word)
baa → class 1 (CP phase word)
```

---

## Dependencies

```
snappy >= 3.3.2   # install via: pip install snappy
latticefit >= 0.2.0
numpy, scipy
sage (for ManifoldHP + exact arithmetic)  # conda install sage
```

---

## Related

- **Website**: [hyperbolicflavorgeometry.org](https://hyperbolicflavorgeometry.org)
- **Papers repo**: [hyperbolic-flavor-geometry](https://github.com/drmlgentry/hyperbolic-flavor-geometry)
- **SSRN**: [ssrn.com/author=11170302](https://ssrn.com/author=11170302)
- **LatticeFit**: [pypi.org/project/latticefit](https://pypi.org/project/latticefit/0.2.0/)

---

*© 2026 Marvin L. Gentry · ORCID: 0009-0006-4550-2663 · drmlgentry@protonmail.com*
