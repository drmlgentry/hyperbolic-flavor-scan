# Hyperbolic Flavor Geometry — Scan Repository

Computational scripts for the spectral universality and slope encoding
analysis of Dehn fillings of the cusped hyperbolic 3-manifolds m003 and m006.

**Associated preprint:** Gentry, M.L. (2026). *A Universal Spectral Phase
Transition in Dehn Fillings of m003 and m006: Lucas Geodesics and Slope
Encoding.* SSRN 6670778.
https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6670778

---

## Key Results

### Phase transition
| Manifold | Threshold n* | log(n*)/log(φ) | Lucas? |
|---|---|---|---|
| m003 | 21 | 6.327 | between L₆ and L₇ |
| m006 | **18** | **6.006** | **exactly L₆** |

For n ≥ n*, log(n) appears as a geodesic length in >80% of all fillings.

### Lucas universality
| n | m003 (N=32) | m006 (N=33) | Lucas index |
|---|---|---|---|
| 18 | 65.6% | 81.8% | L₆ |
| 23 | 96.9% | 87.9% | ≈L₆.₅ |
| 29 | 81.2% | 87.9% | L₇ |
| 37 | 90.6% | 100% | ≈L₇.₅ |
| 47 | 100% | 100% | L₈ |

### Slope encoding (m003, p-coordinate only)
| \|p\| | Baseline | p-fill% | Enrichment |
|---|---|---|---|
| 2 | 6.2% | 20.0% | **3.20×** |
| 3 | 6.2% | 14.3% | **2.29×** |
| 5 | 6.2% | 12.5% | **2.00×** |

**q-coordinate: 0.00× enrichment (zero signal).**

### Kernel rank analysis
Gaussian kernel K_ij = exp(-(ℓ_i - ℓ_j)²/2σ²) on actual geodesic lengths:

| Manifold | N geodesics | R_eff at σ=log(φ) | λ₃/λ₄ gap |
|---|---|---|---|
| M_PMNS = m003(-2,3) | 230 | **3.90** | 1.78 |
| M_CKM = m006(-5,2) | 440 | **3.24** | 2.09 |

Null test: shuffled lengths give identical R_eff (p=0.99); uniform random
gives R_eff≈6.5 (p<10⁻⁴). The hyperbolic length distribution is
intrinsically more compressed than random.

---

## Scripts

### Scan scripts
| Script | Description |
|---|---|
| `scan_flush.py` | Full scan, cutoff=4.5, 32+33 fillings, flushed output |
| `scan_fast.py` | Fast scan, cutoff=4.0, progress indicators |
| `scan_fast_cutoff.py` | Fast scan with timing data per slope |
| `scan_fast2.py` | Fast scan v2 with proper deduplication |

### Analysis scripts
| Script | Description |
|---|---|
| `slope_encoding.py` | Definitive slope encoding analysis (p vs q coordinate) |
| `kernel_sweep.py` | R_eff vs σ sweep for PMNS and CKM geodesic spectra |
| `kernel_rank_analysis.py` | Full kernel rank analysis with null tests |
| `plot_kernel_rank.py` | Generate fig_kernel_rank.pdf/png |

### Lucas structure scripts
| Script | Description |
|---|---|
| `lucas_mode.py` | Lucas mode analysis tools |
| `m006_lucas_scan.py` | Lucas scan for m006 fillings |
| `covering_tower.py` | Covering tower computation |
| `verify_covers_corrected.py` | Verified covering tower results |

### Key manifolds (SnapPy)
```python
import snappy
M_PMNS = snappy.OrientableClosedCensus[1]   # vol=0.98137, H1=Z/5
M_CKM  = snappy.OrientableClosedCensus[43]  # vol=2.02885, H1=Z/5
# Or via Dehn filling:
M_PMNS = snappy.Manifold('m003'); M_PMNS.dehn_fill((-2,3))
M_CKM  = snappy.Manifold('m006'); M_CKM.dehn_fill((-5,2))
# Length spectrum (use .real() not .real):
lengths = [float(x.length.real()) for x in M.length_spectrum(4.5)]
```

---

## Figures
- `fig_kernel_rank.pdf/png` — R_eff vs σ and energy fractions for PMNS and CKM

---

## Requirements
```
snappy>=3.3.2
numpy
matplotlib
scipy
sagemath (for sage kernel)
```
Install: `pip install snappy numpy matplotlib scipy`

---

## Reproducibility
All scans use deterministic slope enumeration with deduplication
`(p,q) ↔ (-p,-q)` only. Cutoff=4.5 for definitive results; 4.0 for
fast scans. Tolerance=0.02 for log(n) matching.

Blacklisted slopes (computation time >10 min at cutoff=4.5):
`(-5,-4), (5,-4), (-4,5), (4,-5)` — timing data in `scan_fast_cutoff.py`.

---

## Related preprints (SSRN)
| Paper | SSRN |
|---|---|
| SU paper (this repo) | [6670778](https://papers.ssrn.com/abstract=6670778) |
| Lucas structure | [6660598](https://papers.ssrn.com/abstract=6660598) |
| HFG Conjectures | [6670398](https://papers.ssrn.com/abstract=6670398) |
| GW phi-lattice | [6669600](https://papers.ssrn.com/abstract=6669600) |
| CKM paper | [6583550](https://papers.ssrn.com/abstract=6583550) |
| PMNS paper | [6583553](https://papers.ssrn.com/abstract=6583553) |

---

## Author
Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663
Independent Researcher, Seattle WA
