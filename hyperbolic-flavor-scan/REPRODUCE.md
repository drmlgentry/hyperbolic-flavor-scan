# Reproducing HFG Results

This repository contains all computational scripts supporting the
Hyperbolic Flavor Geometry (HFG) papers. Each paper has **one canonical
reproduction script** in `reproduce/`. A referee needs to run exactly
one command to verify any paper's claims.

## Requirements

```bash
pip install snappy numpy scipy
```

SnapPy requires Python 3.8–3.11. A conda environment works well:

```bash
conda create -n sage python=3.11
conda activate sage
pip install snappy numpy scipy pandas
```

## The Two Manifolds

All results derive from two compact hyperbolic 3-manifolds:

| Role  | SnapPy name   | Census index              | Volume  | H₁   |
|-------|---------------|---------------------------|---------|------|
| PMNS  | m003(−2,3)    | OrientableClosedCensus[1] | 0.98137 | ℤ/5  |
| CKM   | m006(−5,2)    | OrientableClosedCensus[43]| 2.02885 | ℤ/5  |

**Always load by index**, not by name string.

## Paper Reproduction Scripts

### Paper 1: CKM and PMNS Mixing Matrices

```bash
python reproduce/hfg_reproduce.py
```

Reproduces: PMNS Borel fitness 0.005087 (global minimum),
CKM Gaussian fitness 0.016482, Cabibbo angle error 0.19%.
Runtime: ~2 minutes.

### Paper 2: CP Phase (δ = 195.91°)

```bash
python reproduce/cp_reproduce.py
```

Reproduces: δ = π + φ(aaB) + φ(baa) = 195.91° (0.55% from PDG 197.0°).
Runtime: ~10 seconds.

### Paper 3: Twist Angle Spectrum

```bash
python reproduce/twist_reproduce.py
```

Reproduces all 9 numerical claims in the twist angle census paper:
- CKM isospectrality: φ(aaB)=φ(AbA)=φ(AAb)=92.49° (spread=0)
- δ_CKM = 67.65° from 180−φ(aa) on m006
- θ₂₃_CKM = 2.132° from φ(aaabb) on m006
- PMNS solar angle = 33.62° from 180−φ(abbAB) on m006
- θ₁₂_CKM = 12.64° from 180−φ(AAB) on m003
- mb/mc = 3.2910 from |λ(bbbb)|/|λ(bAbA)| on m003
- δ_CP = 195.91° from π+φ(aaB)+φ(baa) on m003
- MZ/MW = 1.1348 from ell(BBBBB)/ell(baa) on m006
- θ₁₃_ν = 8.546° from 180−φ(BaBBBBB) on m006 [length-7]

Runtime: ~3 minutes.

## Twist Angle Convention

All twist angles in Paper 3 use the **positive-branch eigenvalue**:

    φ(γ) = Im(log λ)   where λ is chosen so Im(log λ) ≥ 0

This gives φ ∈ [0°, 180°] for all loxodromic elements.
Implementation: `G.SL2C(word)` from SnapPy's fundamental group.

The CP phase paper (Paper 2) uses `polished_holonomy()` at 150-bit
precision, which can give negative angles. Both conventions are
consistent: φ_positive = 180° − |φ_polished| when the polished
angle is negative near 180°.

## Verified Canonical Values

| Quantity              | Value      | PDG target | Error  | Words        |
|-----------------------|------------|------------|--------|--------------|
| PMNS Borel fitness    | 0.005087   | —          | —      | aa/aaB/baa   |
| CKM Gaussian fitness  | 0.016482   | —          | —      | aaB/AbA/AAb  |
| δ_CP                  | 195.91°    | 197.0°     | 0.55%  | aaB/baa      |
| δ_CKM                 | 67.65°     | 68.0°      | 0.51%  | aa (m006)    |
| θ₁₂_CKM              | 12.64°     | 13.04°     | 0.31°  | AAB (m003)   |
| θ₂₃_CKM              | 2.132°     | 2.38°      | 0.25°  | aaabb (m006) |
| θ₁₂_ν (solar)        | 33.62°     | 33.41°     | 0.63%  | abbAB (m006) |
| θ₁₃_ν                | 8.546°     | 8.54°      | 0.007° | BaBBBBB(m006)|
| mb/mc                 | 3.2910     | 3.2913     | 0.01%  | bbbb/bAbA    |
| MZ/MW                 | 1.1348     | 1.1345     | 0.025% | BBBBB/baa    |
| CKM isospectral φ     | 92.487°    | —          | —      | aaB=AbA=AAb  |

## Repository Structure

```
reproduce/          ← START HERE (one script per paper)
  hfg_reproduce.py  ← CKM + PMNS fitness
  cp_reproduce.py   ← CP phase 195.91°
  twist_reproduce.py← twist angle census claims

archive/            ← development scans (not needed for reproduction)
  scans/            ← all intermediate scan scripts
  data/             ← CSV output files from scans

README.md           ← overview
REPRODUCE.md        ← this file
```

## SSRN Preprints

| SSRN     | Paper                          |
|----------|--------------------------------|
| 6775158  | Unified HFG (all results)      |
| 6583550  | CKM matrix                     |
| 6583549  | PMNS matrix                    |
| 6583551  | CP phase                       |
| 6583553  | Twist angle spectrum           |

## Contact

Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663
