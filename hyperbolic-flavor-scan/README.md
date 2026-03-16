# hyperbolic-flavor-scan

Numerical scan code and census data for the hyperbolic flavor geometry program.
Companion to [hyperbolic-flavor-geometry](https://github.com/drmlgentry/hyperbolic-flavor-geometry).

## Manifold identifications

| SnapPy name | Census index | Mathematical name | Volume | H1 |
|---|---|---|---|---|
| m003 | OrientableClosedCensus[1] | Meyerhoff manifold | 0.9814 | Z/5 |
| m006 | OrientableClosedCensus[43] | — | 2.0289 | Z/5 |

Note: OrientableClosedCensus[0] is the **Weeks manifold** (vol=0.9427, H1=(Z/5)^2),
a distinct manifold. Published lambda1 values (Inoue 2001: ~27.8) apply to the
Weeks manifold, not to m003 or m006.

## Directory structure
```
scans/      Canonical scan scripts
analysis/   Post-processing, table generation, figures
figures/    Output PDFs and PNGs
data/       CSV result files
archive/    Development/debug scripts (not canonical)
```

## Key scripts

| Script | Purpose | Paper |
|---|---|---|
| `scans/word_triple_scan_corrected.py` | CKM word triple scan on m006 | Paper 5 (PRD) |
| `scans/pmns_borel_scan_v2.py` | PMNS Borel construction on m003 | Paper 6 (PRD) |
| `scans/twist_census.py` | A-factor twist census, lengths 1-5 | Paper 8 (NPB) |
| `scans/selberg_setup.py` | Selberg zeta zero search, lambda1 estimates | Paper 8 (NPB) |
| `scans/check_lambda1_v2.py` | Manifold identification + eigenvalue context | Paper 8 (NPB) |
| `analysis/paper8_analysis.py` | Census deduplication + SM coincidence tables | Paper 8 (NPB) |
| `scans/verify_all.py` | Reproduce all key results | All papers |

## Key results

### CKM (m006)
- Best triple: aaB/AbA/AAb | fitness F=0.01729 | J≈0
- phi(b)=89.16 deg ≈ pi/2 explains triple isospectrality and J_CKM≈0
- 180-phi(aa) = 67.65 deg ~ delta_CKM = 68.0 deg (0.5% error)

### PMNS (m003 = Meyerhoff manifold)
- Best triple: aa/ab/aB | column perm (0,2,1) | fitness F=0.01897
- phi_aa - phi_ab + phi_aB = 203.5 deg vs PDG delta_CP=197 deg (3.3%, 0 free params)
- |lam(bbbb)|/|lam(bAbA)| = 3.2910 ~ m_b/m_c = 3.2913 (0.01%, MSbar scheme)

### Selberg eigenvalue estimates (Paper 8, first for these manifolds)
- lambda1(m003/Meyerhoff) ~ 2.48 (leading-term estimate)
- lambda1(m006) ~ 2.82 (leading-term estimate)
- Both >> 1/4; consistent with Selberg eigenvalue conjecture

## Quickstart
```bash
pip install snappy scipy numpy pandas matplotlib --break-system-packages
python scans/verify_all.py
```

## Author

Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663
