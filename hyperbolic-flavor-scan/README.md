# Hyperbolic Flavor Geometry

Computational framework for the **Hyperbolic Flavor Geometry (HFG)** programme — deriving Standard Model flavor parameters (CKM/PMNS mixing matrices, CP violation, fermion masses) from the arithmetic geometry of compact hyperbolic 3-manifolds.

**Website:** [hyperbolicflavorgeometry.org](https://hyperbolicflavorgeometry.org)

## Key Results

| Result | Manifold | Accuracy |
|---|---|---|
| PMNS lepton mixing | m003(-2,3) = Meyerhoff manifold | fitness 0.005087 |
| CKM quark mixing | m006(-5,2) | fitness 0.016482 |
| CP phase delta = 195.91 deg | holonomy of m003 | 0.55% vs PDG 197.0 |
| m_mu/m_e = 208 | Eisenstein norm N(16+12*omega) | 0.59% |
| m_tau/m_e = 3477 | Eisenstein norm N(68+37*omega) | 0.006% |
| Gal(m003) = Z/2 = Weyl(SU(2)) | cusp field x^2-x+1 | exact |
| Gal(m006) = S3 = Weyl(SU(3)) | cusp field x^3+2x+1 | exact |
| Gal(m019) = S4 = Weyl(SU(4)) | cusp field x^4-x-1 | exact |
| Dual surgery: m003(-2,3) = m019(2,1) | Meyerhoff manifold | 15 sig. figs. |
| Delta invariant delta(m019)=12, delta(m178)=34 | disc=-283 cusp field | exact |
| Lucas trace: 2*cosh(2m*log(phi)) = L_{2m} | golden ratio identity | exact |
| vol(m019) = 3*vol(M_PMNS) | Bloch-Wigner: D(z)=v0 all tetrahedra | 2e-12 |
| vol(m178) = 4*vol(M_PMNS) | PSL(2,Z) orbit period 3 | 2e-12 |

## Canonical Manifolds

- **m003** — SU(2) cusped parent, cusp shape e^{i*pi/3}, trace field Q(sqrt(-3))
- **m006** — SU(3) cusped parent, cusp shape satisfies x^3+2x+1, disc=-59
- **m019** — SU(4) cusped parent, cusp shape satisfies x^4-x-1, disc=-283
- **m003(-2,3) = m019(2,1)** — Meyerhoff manifold M_PMNS (minimum-volume closed hyperbolic 3-manifold), vol = 0.9813688289...
- **m006(-5,2)** — M_CKM
- **m178** — disc=-283 sibling of m019; vol = 4*vol(M_PMNS)

## Bloch Volume Quantum (New 2026-06-01)

Every tetrahedral shape parameter in m019 and m178 satisfies D(z) = vol(M_PMNS), where D is the Bloch-Wigner dilogarithm. The three distinct shapes form a single PSL(2,Z) orbit of period 3 under T(z) = 1/(1-z):

```
z_A -> z_B -> z_S -> z_A   (errors < 2e-15)
D(z_A) = D(z_B) = D(z_S) = vol(M_PMNS) = 0.9813688289...
vol(m019) = 3 * vol(M_PMNS)
vol(m178) = 4 * vol(M_PMNS)
```

## Reproduction

Requires [SnapPy](http://snappy.computop.org) and [SageMath](https://www.sagemath.org).

```python
import snappy

# Verify dual surgery
M1 = snappy.Manifold("m003(-2,3)")
M2 = snappy.Manifold("m019(2,1)")
print(M1.volume())              # 0.9813688289...
print(M2.volume())              # 0.9813688289...
print(M1.is_isometric_to(M2))   # True

# Verify Bloch volume quantum
M19  = snappy.Manifold("m019")
M178 = snappy.Manifold("m178")
v0   = float(M1.volume())
print(float(M19.volume())  / v0)   # 3.0000000000
print(float(M178.volume()) / v0)   # 4.0000000000
```

**SnapPy note:** Always use `OrientableClosedCensus[1]` (M_PMNS) and `OrientableClosedCensus[43]` (M_CKM) by index.

## Papers

| SSRN | Title | Journal |
|---|---|---|
| [6859979](https://ssrn.com/abstract=6859979) | Bloch Volume Quantum | Experimental Mathematics (submitted) |
| [6854378](https://ssrn.com/abstract=6854378) | Lucas Trace Identity | MRL (under review) |
| [6851440](https://ssrn.com/abstract=6851440) | Delta Invariant, disc=-283 | AGT (under review) |
| [6848478](https://ssrn.com/abstract=6848478) | Mu-Function | — |
| [6845778](https://ssrn.com/abstract=6845778) | Dual Surgery / Pati-Salam | Proc. AMS (under review) |
| [6840418](https://ssrn.com/abstract=6840418) | BPS Lepton Masses on X0(11) | LMP (under review) |
| [6840324](https://ssrn.com/abstract=6840324) | Gauss Polynomial / WRT | — |
| [6840322](https://ssrn.com/abstract=6840322) | Galois-Gauge Correspondence | CNTP (submitted) |
| [6815721](https://ssrn.com/abstract=6815721) | Level-11 Automorphic Structure | — |
| [6775158](https://ssrn.com/abstract=6775158) | Unified HFG Framework | — |

## Python Package

```bash
pip install latticefit
```

[LatticeFit v0.3.1 on PyPI](https://pypi.org/project/latticefit/) — scan hyperbolic manifold holonomy for flavor parameter fitness.

## Author

Marvin L. Gentry — Independent Researcher, Seattle WA
ORCID: [0009-0006-4550-2663](https://orcid.org/0009-0006-4550-2663)
Email: drmlgentry@protonmail.com
Website: [hyperbolicflavorgeometry.org](https://hyperbolicflavorgeometry.org)
