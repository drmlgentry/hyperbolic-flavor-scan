# Hyperbolic Flavor Geometry — Papers

**Author:** Marvin L. Gentry  
**Email:** drmlgentry@protonmail.com  
**ORCID:** 0009-0006-4550-2663  
**Last updated:** May 17, 2026

## Active Submission Portfolio

| ID | Journal | Title | Status |
|---|---|---|---|
| DS14327 | PRD | Unified HFG | With editors |
| PLB-D-26-01341 | PLB | CP phase (twist angles) | With editor |
| PLB-D-26-01006 | PLB | Charge conjugation | Under review |
| AHPO-D-26-00255 | AHP | PMNS mixing | New submission |
| AHPO-D-26-00231 | AHP | CP holonomy | New submission |
| JGP13076 | JGP | CKM mixing | Submitted |
| TRGR-D-26-00059 | Trieste | Discrete mixing | Editor assigned |
| JMP26-AR-01272 | JMP | Qubit gates | Active |

## Final Paper Files (May 17, 2026)

### Core HFG Papers

| File | Description | Submitted to |
|---|---|---|
| `hfg-unified/gentry-hfg-unified.tex` | Unified HFG (all 4 pillars) | PRD DS14327 |
| `hyperbolic-flavor-ckm-v2/gentry-ckm-final.tex` | CKM matrix from m006 | JGP JGP13076 |
| `hyperbolic-flavor-pmns/gentry-pmns-final.tex` | PMNS matrix from m003 | AHP AHPO-D-26-00255 |
| `hyperbolic-flavor-cp-v2/gentry-cp-final.tex` | CP phase from twist angles | PLB PLB-D-26-01341 |

### Figures

| File | Description |
|---|---|
| `hyperbolic-flavor-ckm-v2/ckm_combined_figure.pdf` | CKM sphere + triangle |

## Canonical Numerical Values

All values verified against `hfg_reproduce.py` ground truth output.

**PMNS (m003(-2,3), OrientableClosedCensus[1]):**
- Borel fitness: **0.005087** (global minimum, 15 word triples all agree)
- CP phase: **195.91°** (PDG 197.0°, 0.55% error, zero free parameters)
- Trace field: **ℚ(√−3)** imaginary quadratic
- Lepton norms: μ→208 (0.59%), τ→3477 (0.006%) in ℤ[ω]
- Tower: {2,3,11} = {L₀,L₂,L₅} through degree 19

**CKM (m006(-5,2), OrientableClosedCensus[43]):**
- Gaussian fitness: **0.016482** (σ=0.49 fixed)
- Cabibbo angle error: **0.19%**
- Trace field: **ℚ(√17)** real quadratic, tr(ρ(aa))=3−√17
- Quark norms: all 6 in ℤ[√17], p<0.002
- Tower: {11} = {L₅} through degree 19

**Note on covering tower:** The earlier claim {2,3,7,11,29} was a
conflation with the geodesic length spectrum (SU paper). Verified
covering tower torsion primes through degree 19: {2,3,11} (PMNS)
and {11} (CKM). Primes 7 and 29 do not appear through degree 19.

## Reproducibility

All results reproducible via:
```bash
conda run -n sage python hfg_reproduce.py
```

Full verification suite in `hyperbolic-flavor-scan` repository.

## SSRN Preprints

| SSRN ID | Title |
|---|---|
| 6775158 | Unified HFG (submitted May 16 2026) |
| 6583550 | CKM mixing |
| 6583549 | PMNS mixing |
| 6583553 | CP A-factor |
| 6754501 | Lucas structure (needs correction: tower primes) |
| 6761981 | Covering tower |

## arXiv

Endorsement code obtained (SIQ). Endorsement request sent to Kofman
(May 2026). Second endorsement request pending.
