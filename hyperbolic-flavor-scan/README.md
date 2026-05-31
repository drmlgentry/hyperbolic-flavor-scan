# Hyperbolic Flavor Geometry

Computational framework for the **Hyperbolic Flavour Geometry (HFG)** programme — deriving Standard Model flavor parameters (CKM/PMNS mixing matrices, CP violation, fermion masses) from the arithmetic geometry of compact hyperbolic 3-manifolds.

## Key Results

| Result | Manifold | Accuracy |
|---|---|---|
| PMNS lepton mixing | m003(-2,3) = Meyerhoff manifold | fitness 0.005087 |
| CKM quark mixing | m006(-5,2) | fitness 0.016482 |
| CP phase δ = 195.91° | holonomy of m003 | 0.55% vs PDG 197.0° |
| m_μ/m_e = 208 | Eisenstein norm N(16+12ω) | 0.59% |
| m_τ/m_e = 3477 | Eisenstein norm N(68+37ω) | 0.006% |
| Galois–Weyl: Gal(m003)=Z/2=Weyl(SU(2)) | cusp field x²-x+1 | exact |
| Galois–Weyl: Gal(m006)=S3=Weyl(SU(3)) | cusp field x³+2x+1 | exact |
| Galois–Weyl: Gal(m019)=S4=Weyl(SU(4)) | cusp field x⁴-x-1 | exact |
| Dual surgery: m003(-2,3) = m019(2,1) | Meyerhoff manifold | 15 sig. figs. |
| Delta invariant δ(m019)=12, δ(m178)=34 | disc=-283 cusp field | exact |
| Lucas trace: 2cosh(2m·log φ) = L_{2m} | golden ratio identity | exact |

## Canonical Manifolds

- **m003** — SU(2) cusped parent, trace field Q(√-3), cusp shape e^{iπ/3}
- **m006** — SU(3) cusped parent, trace field Q(√-59), cusp shape satisfies x³+2x+1
- **m019** — SU(4) cusped parent, trace field Q(√-283), cusp shape satisfies x⁴-x-1
- **m003(-2,3)** = **m019(2,1)** — Meyerhoff manifold M_PMNS (minimum-volume closed hyperbolic 3-manifold)
- **m006(-5,2)** — M_CKM

## Reproduction

Requires [SnapPy](http://snappy.computop.org) and [SageMath](https://www.sagemath.org).

```python
# Verify dual surgery
import snappy
M1 = snappy.Manifold("m003(-2,3)")
M2 = snappy.Manifold("m019(2,1)")
print(M1.volume())   # 0.9813688289...
print(M2.volume())   # 0.9813688289...
print(M1.is_isometric_to(M2))  # True
```

Canonical reproduction script: `hyperbolic-flavor-scan/hfg_reproduce.py`

**Note:** Always use `OrientableClosedCensus[1]` (PMNS) and `OrientableClosedCensus[43]` (CKM) by index in SnapPy. Do not load by name string "m006" (two distinct census manifolds share that name).

## Papers (preprints on SSRN)

| SSRN | Title |
|---|---|
| [6775158](https://ssrn.com/abstract=6775158) | Unified HFG framework |
| [6840322](https://ssrn.com/abstract=6840322) | Galois–Gauge correspondence (SU(2), SU(3), SU(4)) |
| [6840418](https://ssrn.com/abstract=6840418) | Lepton masses as BPS states on X₀(11) |
| [6845778](https://ssrn.com/abstract=6845778) | Dual surgery / Pati-Salam (Proc. AMS submitted) |
| [6851440](https://ssrn.com/abstract=6851440) | Delta invariant for disc=-283 manifolds (AGT submitted) |
| [6854378](https://ssrn.com/abstract=6854378) | Lucas trace identity (MRL submitted) |

## Python Package

[![PyPI](https://img.shields.io/pypi/v/latticefit)](https://pypi.org/project/latticefit/)

```bash
pip install latticefit
```

LatticeFit v0.2.0: scan hyperbolic manifold holonomy for flavor parameter fitness.

## Author

Marvin L. Gentry — Independent Researcher, Seattle WA  
ORCID: [0009-0006-4550-2663](https://orcid.org/0009-0006-4550-2663)  
Email: drmlgentry@protonmail.com
