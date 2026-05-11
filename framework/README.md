# Hyperbolic Flavor Geometry — Framework Repository

LaTeX sources and code for the Hyperbolic Flavor Geometry (HFG) research
program: a geometric explanation of Standard Model flavor structure using
compact hyperbolic 3-manifolds.

**Program overview:** The HFG framework derives fermion masses, CKM/PMNS
mixing matrices, and CP violation from the hyperbolic geometry of two
closed manifolds M_PMNS = m003(-2,3) and M_CKM = m006(-5,2), via
holonomy representations and Gaussian kernel overlap matrices on S².

---

## Paper Portfolio

### Submitted to journals (15 submissions, 8 journals)

| Paper | Journal | ID | SSRN |
|---|---|---|---|
| Chirality | PLB | PLB-D-26-01006 | — |
| Neutrino | JPhysG | JPhysG-105833 | [6631218](https://papers.ssrn.com/abstract=6631218) |
| Lucas structure | JNT | JNTH-D-26-00538 | [6660598](https://papers.ssrn.com/abstract=6660598) |
| Covering tower | MMJ | submitted | — |
| Weeks/Dehn | LMP | MATH-D-26-00263 | — |
| Mixing Operators | LMP | MATH-D-26-00165 | — |
| CP Phases | LMP | MATH-D-26-00166 | — |
| CKM | RiP | RINP-D-26-00327 | [6583550](https://papers.ssrn.com/abstract=6583550) |
| PMNS | RiP | RINP-D-26-00328 | [6583553](https://papers.ssrn.com/abstract=6583553) |
| CP A-factor | RiP | RINP-D-26-00329 | [6583551](https://papers.ssrn.com/abstract=6583551) |
| Twist Spectrum | RiP | RINP-D-26-00330 | [6583549](https://papers.ssrn.com/abstract=6583549) |
| Homology Asymmetry | AIF | 2026120 | — |
| Alexander polynomial | JKTR | JKTR-S-26-00044 | — |

### Preprints (SSRN only)

| Paper | SSRN | Posted |
|---|---|---|
| SU paper: Spectral Universality | [6670778](https://papers.ssrn.com/abstract=6670778) | Apr 28, 2026 |
| HFG Conjectures (dark sector, cosmology) | [6670398](https://papers.ssrn.com/abstract=6670398) | Apr 28, 2026 |
| GW phi-lattice / dark matter | [6669600](https://papers.ssrn.com/abstract=6669600) | Apr 28, 2026 |

---

## Repository Structure

```
framework/
├── papers/
│   ├── spectral-universality/   SU paper (long + short, tex + pdf)
│   ├── lucas-structure/         Lucas bridge theorem (JNT submitted)
│   ├── hfg-conjectures/         Conjectures on dark sector + cosmology
│   ├── hyperbolic-flavor-ckm/   CKM paper (RiP submitted)
│   ├── hyperbolic-flavor-pmns/  PMNS paper (RiP submitted)
│   ├── hyperbolic-flavor-cp/    CP paper (EPJC/RiP)
│   ├── hyperbolic-flavor-twist/ Twist spectrum paper
│   ├── covering-tower/          Covering tower paper (MMJ)
│   ├── holonomy-cp/             Holonomy CP paper
│   ├── neutrino/                Neutrino paper (JPhysG)
│   ├── chirality/               Chirality paper (PLB)
│   ├── weeks-dehn/              Weeks/Dehn paper (LMP)
│   └── lucas-structure/         Lucas structure (JNT)
├── figures/                     Shared figures directory
│   ├── fig_effective_rank.pdf/png
│   ├── fig_kernel_rank.pdf/png
│   ├── fig_two_tier_schematic.pdf/png
│   ├── fig_mass_lattice.png
│   └── fig_phase_transition.png
└── code/
    ├── verify_ckm.py
    ├── verify_mixing.py
    └── mass_only_eval.py
```

---

## Core Mathematical Results

### The Lucas-geodesic bridge (exact)
ℓ = k·log(φ) ⟺ |tr(γ)| = Lₖ = φᵏ + φ⁻ᵏ

### Three-generation threshold (proved algebraically)
1/6 < (log φ)² < 1/3

Consequently, σ = log(φ) supports exactly 3 active spherical harmonic
modes on S² — a geometric realisation of three fermion generations.

### Phase transition (computational)
- n*(m006) = 18 = L₆ exactly (log 18/log φ = 6.006)
- n*(m003) = 21 (between L₆ and L₇)

### Slope encoding
p-coordinate encodes log|p| geodesics: 3.20× (|p|=2), 2.29× (|p|=3),
2.00× (|p|=5). q-coordinate: zero enrichment.

### Kernel rank at σ = log(φ)
- M_PMNS: R_eff = 3.90
- M_CKM:  R_eff = 3.24

---

## Key Constants
```python
phi = (1 + sqrt(5)) / 2        # = 1.6180339...
log_phi = log(phi)              # = 0.4812118...  fundamental geodesic unit
sigma_opt = (3/2) * log_phi    # = 0.7218178...  PMNS fitting parameter
# L_k = phi^k + phi^(-k):
# L_0=2, L_1=1, L_2=3, L_3=4, L_4=7, L_5=11, L_6=18, L_7=29, L_8=47
```

---

## Related Repositories
- [hyperbolic-flavor-scan](https://github.com/drmlgentry/hyperbolic-flavor-scan)
  — SnapPy scan scripts and geodesic analysis
- [latticefit](https://pypi.org/project/latticefit/)
  — PyPI package for φ-lattice mass fitting (v0.2.0, Lucas mode)

---

## Author
Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663
Independent Researcher, Seattle WA
arXiv endorsement pending: SIQW6F (math.GT)


## Paper Status (May 2026)

| Paper | Journal | ID | Status |
|-------|---------|-----|--------|
| CKM (corrected) | Results in Physics | RINP-D-26-00327 | Correction submitted |
| PMNS (corrected) | Results in Physics | RINP-D-26-00328 | Correction submitted |
| CP Phases | Results in Physics | RINP-D-26-00329 | Under review |
| Twist Spectrum | Results in Physics | RINP-D-26-00330 | Under review |
| Weeks/Dehn | JGP | JGP13023 | With editor |
| Mixing Operators | Transformation Groups | TRGR-D-26-00059 | Submitted |
| CP Phases (math) | Annales Henri Poincare | AHPO-D-26-00231 | Submitted |
| Qubit Gates | JMP | JMP26-AR-01272 | With editor |
| Chirality | PLB | PLB-D-26-01006 | Under review |
| Neutrino | J. Phys. G | JPhysG-105833 | Submitted |
| Alexander poly | JKTR | JKTR-S-26-00044 | With editor |
| Covering Tower | MMJ | -- | Submitted |
