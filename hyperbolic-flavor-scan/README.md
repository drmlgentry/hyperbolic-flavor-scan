# hyperbolic-flavor-scan

Numerical scan code and census data for the hyperbolic flavor geometry program.
Companion repository to [hyperbolic-flavor-geometry](https://github.com/drmlgentry/hyperbolic-flavor-geometry).

## What this repo contains

Scan scripts, analysis tools, and result data for all eight papers in the program.
The central result: compact hyperbolic 3-manifolds **m003** and **m006** (both with H1=Z/5)
reproduce the CKM and PMNS mixing matrices and encode the full SM flavor parameter
set in their A-factor twist angle spectra.

## Directory structure
```
scans/          Canonical scan scripts (use these for reproduction)
analysis/       Post-processing, table generation, figure scripts
figures/        Output figures (PDF + PNG)
data/           CSV result files from completed scans
archive/        Development/debug scripts (not canonical)
```

## Key scripts

| Script | Purpose | Paper |
|---|---|---|
| `scans/word_triple_scan_corrected.py` | CKM word triple scan on m006 | Paper 5 (PRD) |
| `scans/pmns_borel_scan_v2.py` | PMNS Borel construction scan on m003 | Paper 6 (PRD) |
| `scans/twist_census.py` | Full A-factor twist census (m003+m006, length 1-5) | Paper 8 (NPB) |
| `scans/selberg_setup.py` | Selberg zeta zero search, lambda1 estimates | Paper 8 (NPB) |
| `analysis/paper8_analysis.py` | Census deduplication and SM coincidence table | Paper 8 (NPB) |
| `analysis/results_summary.py` | Side-by-side CKM+PMNS fitness summary | Papers 5-6 |

## Key results

### CKM (m006, OrientableClosedCensus[43])
- Best triple: aaB/AbA/AAb | fitness F=0.01729 | J≈0
- Generator twist: φ(b)=89.16° ≈ π/2 (explains triple isospectrality)

### PMNS (m003, OrientableClosedCensus[1])
- Best triple: aa/ab/aB | column perm (0,2,1) | fitness F=0.01897
- CP phase: φ_aa − φ_ab + φ_aB = 203.5° vs PDG 197° (3.3% error, 0 free params)

### Twist angle spectrum (Paper 8)
- m006 aa: 180°−φ = 67.65° ~ δ_CKM = 68.0° (0.51% error)
- m006 abbAB: 180°−φ = 33.62° ~ θ12_ν = 33.41° (0.63% error)
- m003 |λ(bbbb)|/|λ(bAbA)| = 3.2910 ~ m_b/m_c = 3.2913 (0.01% error, MS-bar scheme)
- Selberg estimates: λ1(m003)≈2.48, λ1(m006)≈2.82

## Quickstart
```bash
pip install snappy scipy numpy pandas matplotlib --break-system-packages
python scans/verify_all.py        # reproduce all key results
python scans/twist_census.py      # regenerate twist census CSVs
python analysis/paper8_analysis.py  # regenerate Paper 8 coincidence tables
```

## Citation

If you use this code, please cite the relevant papers from the
[hyperbolic-flavor-geometry](https://github.com/drmlgentry/hyperbolic-flavor-geometry) repository.

## Author

Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663
