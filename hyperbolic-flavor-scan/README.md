# hyperbolic-flavor-scan

Numerical scan code and census data for the hyperbolic flavor geometry program.
Companion to [hyperbolic-flavor-geometry](https://github.com/drmlgentry/hyperbolic-flavor-geometry).

## Manifold identifications

| SnapPy name | Census index | Mathematical name | Volume | H1 |
|---|---|---|---|---|
| m003 | OrientableClosedCensus[1] | Meyerhoff manifold | 0.9814 | Z/5 |
| m006 | OrientableClosedCensus[43] | — | 2.0289 | Z/5 |

**Note:** OrientableClosedCensus[0] is the Weeks manifold (vol=0.9427, H1=(Z/5)^2),
a distinct manifold. Published λ1 values (Inoue 2001: ~27.8) apply to the Weeks
manifold, not to m003 or m006. Our λ1 estimates are the first for these manifolds.

## Directory structure
```
scans/      Canonical scan scripts and verification entry point
analysis/   Figure generation scripts
figures/    Output PDFs and PNGs (all papers)
data/       CSV result files from completed scans
archive/    Development/debug scripts (not canonical)
```

## Key scripts

| Script | Purpose | Paper |
|---|---|---|
| `scans/word_triple_scan_corrected.py` | CKM word triple scan on m006 | Paper 5 (PRD) |
| `scans/pmns_borel_scan_v2.py` | PMNS Borel N-factor construction on m003 | Paper 6 (PRD) |
| `scans/twist_census.py` | A-factor twist census, lengths 1-5 | Paper 8 |
| `scans/twist_census_len7.py` | Extended census to length 7 | Paper 8 |
| `scans/selberg_setup.py` | Selberg zeta zero search, λ1 estimates | Paper 8 |
| `scans/paper8_analysis.py` | Census deduplication + SM coincidence tables | Paper 8 |
| `scans/spectrum_analysis.py` | Spectral gap analysis (length-7 result) | Paper 8 |
| `scans/verify_all.py` | Reproduce all key results | All papers |

## Key results

### CKM (m006)
- Best triple: aaB/AbA/AAb | fitness F=0.01729 | J≈0
- φ(b)=89.16° explains triple isospectrality and J_CKM≈0
- 180°−φ(aa)=67.65° ~ δ_CKM=68.0° (0.51% error)

### PMNS (m003 = Meyerhoff manifold)
- Best triple: aa/ab/aB | column perm (0,2,1) | fitness F=0.01897
- φ_aa − φ_ab + φ_aB = 203.5° vs PDG δ_CP=197° (3.3%, 0 free params)
- |λ(bbbb)|/|λ(bAbA)| = 3.2910 ~ m̄_b/m̄_c = 3.2913 (0.01%, MS-bar scheme)

### Twist angle spectrum (Papers 7-8)
- m006 word abbAB: 180°−φ = 33.62° ~ θ12_ν=33.41° (0.63%)
- m006 length-6 near-identity geodesic: φ_fold = 0.0053° (word aaabab)
- Spectral gap [0.005°, 1.61°] on m006 at lengths ≤ 7; θ13_CKM=0.201° falls inside

### Selberg eigenvalue estimates (first for these manifolds)
- λ1(m003/Meyerhoff) ≈ 2.48 (leading-term estimate)
- λ1(m006) ≈ 2.82 (leading-term estimate)

## Figures

| File | Description | Paper |
|---|---|---|
| `figures/fig_poincare_ball.pdf` | Geodesic axes in Poincaré ball | Paper 5 (CKM) |
| `figures/fig_spherical_triangle.pdf` | Spherical triangle of axis angles | Paper 5 (CKM) |
| `figures/sigma_sensitivity.png` | Fitness vs σ parameter | Paper 5 (CKM) |
| `figures/fig_cp_pareto.png` | Pareto frontier F vs J | Paper 7 (CP) |
| `figures/fig_cp_afactor.png` | A-factor bar chart | Paper 7 (CP) |
| `figures/fig_twist_spectrum.pdf` | Twist angle spectrum φ vs word length | Paper 8 |

## Quickstart
```bash
pip install snappy scipy numpy pandas matplotlib --break-system-packages
python scans/verify_all.py
```

## Author

Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663
