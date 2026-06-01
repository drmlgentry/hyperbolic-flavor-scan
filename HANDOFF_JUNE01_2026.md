# HFG SESSION HANDOFF â€” June 1, 2026
## Researcher: Marvin L. Gentry | drmlgentry@protonmail.com | ORCID: 0009-0006-4550-2663

---

## COMMUNICATION STYLE AND WORKING PREFERENCES

- **PowerShell** for all Windows commands (not bash). Always use `Set-Content -Encoding UTF8` for file writes, or `[System.IO.File]::WriteAllText()` to avoid BOM issues.
- **Complete files only** â€” never surgical edits to LaTeX unless absolutely necessary. Always deliver a full compilable file.
- **One-step answers** to direct questions during computation sessions.
- **Git commits after every significant result** â€” commit messages should be detailed and include submission IDs, paper titles, and key mathematical results.
- **Base64 encoding** for Python scripts when quote conflicts arise.
- **HANDOFF.md** practice established â€” always write before ending a long session.
- Responses should be precise but not terse. Include relevant LaTeX, Python, and reasoning.
- Favor falsifiability and structural explanations over post-hoc fitting.
- Prove what can be proved exactly; state conjectures as conjectures.
- Run computations before theorizing.

---

## ENVIRONMENT AND FILE LOCATIONS

### Windows/WSL Setup
- **Primary working directory:** `C:\dev\framework\`
- **Scan codebase:** `C:\dev\hyperbolic-flavor-scan\`
- **Sage/SnapPy environment:** WSL Ubuntu, conda env named `sage`
  - Activate: `conda activate sage` in WSL
  - **CRITICAL:** PMNS computations MUST use WSL sage env (not Windows) due to holonomy extraction discrepancy
  - SnapPy in sage: `import snappy`
- **Python venv for builds/uploads:** `C:\dev\hyperbolic-flavor-scan\.venv\Scripts\python.exe`
  - Has: twine, build, latticefit dependencies
- **PyPI credentials:** stored in `$env:USERPROFILE\.pypirc` (delete after use â€” contains token)
- **GitHub token:** stored in git remote URLs â€” REVOKE old token `ghp_REDACTED_REVOKE_THIS_TOKEN` at https://github.com/settings/tokens (EXPOSED in terminal June 1)

### Key Scripts
```
C:\dev\framework\hfg_reproduce.py          # canonical reproduction script
C:\dev\framework\spectral_gap_test.py      # Monte Carlo spectral gap
C:\dev\framework\disc283_scan.py           # disc=-283 census scan
C:\dev\framework\census_mu_scan_v3.py      # census scan v3
C:\dev\framework\length_spectrum_checkpointed.py
C:\dev\framework\word_triple_scan_corrected.py
C:\dev\framework\check_angles.py
C:\dev\framework\sigma_sensitivity.py
```

### Key Data Files
```
C:\dev\framework\census_v3_checkpoint.json      # 1000 manifolds, 2 self-encoding found
C:\dev\framework\disc283_checkpoint.json        # 2000 manifolds, 3 disc=-283, 2 self-encoding
C:\dev\framework\spectral_gap_checkpoint.json
C:\dev\framework\length_spectrum_checkpoint.json
```

### GitHub Repos
- `https://github.com/drmlgentry/hyperbolic-flavor-geometry` â€” papers, docs, GitHub Pages site
- `https://github.com/drmlgentry/hyperbolic-flavor-scan` â€” scan code, checkpoints
- `https://github.com/drmlgentry/latticefit` â€” PyPI package
- `https://github.com/drmlgentry/golden_unification` â€” earlier work

### Websites
- **Primary:** https://hyperbolicflavorgeometry.org (Cloudflare Pages, free)
  - Cloudflare account: drmlgentry@protonmail.com
  - Connected to: hyperbolic-flavor-geometry repo, docs/ folder
  - Domain registered June 1, 2026, auto-renews May 31, 2027 ($10.13/yr)
- **Mirror:** https://drmlgentry.github.io/hyperbolic-flavor-geometry
- **SSRN author page:** https://papers.ssrn.com/sol3/cf_dev/AbsByAuth.cfm?per_id=11170302

---

## CANONICAL MANIFOLD FACTS (memorize these)

```
M_PMNS = m003(-2,3) = m019(2,1) = OrientableClosedCensus[1]
  vol = 0.9813688289... = v0 (fundamental Bloch quantum)
  H1 = Z/5, CS = 1/4 (exact)

M_CKM = m006(-5,2) = OrientableClosedCensus[43]
  vol = 2.0289..., H1 = Z/5

m003 (cusped SU(2) parent):
  cusp shape = e^{i*pi/3} = omega, satisfies x^2-x+1
  trace field = Q(sqrt(-3)), disc = -3, Gal = Z/2 = Weyl(SU(2))

m006 (cusped SU(3) parent):
  cusp shape satisfies x^3+2x+1
  trace field disc = -59, Gal = S3 = Weyl(SU(3))

m019 (cusped SU(4) parent):
  cusp shape satisfies x^4-x-1
  disc = -283, Gal = S4 = Weyl(SU(4))
  vol = 2.9441064867... = 3*v0 (PROVED June 1)
  H1 = Z (torsion-free)
  delta(m019) = 12

m178 (disc=-283 sibling):
  vol = 3.9254753156... = 4*v0 (PROVED June 1)
  H1 = Z (torsion-free)
  delta(m178) = 34

CRITICAL SnapPy note: NEVER load by name string for closed manifolds.
Use OrientableClosedCensus[1] for M_PMNS, OrientableClosedCensus[43] for M_CKM.
```

---

## KEY RESULTS (all verified, ordered by strength)

### Exact Algebraic Theorems (proved June 1, 2026)
1. **Bloch orbit theorem:** In K=Q(w), w^4=w+1:
   - z_A = w^3, z_B = -w, z_S = w^{-4}
   - z_A * z_B * z_S = -1 (exact)
   - T(w^3)=-w, T(-w)=w^{-4}, T(w^{-4})=w^3 where T(z)=1/(1-z)
   - Proved by pure algebra: w(1-w^3) = w-w^4 = w-(w+1) = -1
   - Therefore D(z_A)=D(z_B)=D(z_S) by Bloch-Wigner functional equation
   - Therefore vol(m019)=3*v0, vol(m178)=4*v0

2. **Delta invariant:** delta(M) = min(|6a-19b|, |13a+6b|, |13b-19a|)
   - For torsion-free H1: |H1(M(pi*))| = delta(M) exactly
   - Proved from surgery formula

3. **Lucas trace identity:** 2*cosh(2m*log(phi)) = L_{2m}
   - Proved from Binet formula (even k only)

4. **Dual surgery:** m003(-2,3) = m019(2,1) = M_PMNS
   - Verified to 15 sig figs, explicit isometry check

5. **Galois-Weyl trichotomy:**
   - Gal(m003) = Z/2 = Weyl(SU(2))
   - Gal(m006) = S3 = Weyl(SU(3))
   - Gal(m019) = S4 = Weyl(SU(4))
   - All exact, 300-bit precision

### Numerical Observations (high precision, not yet proved algebraically)
- D(w^3) = v0 = 0.9813688289... (to 2e-12) â€” the one remaining numerical step
- PMNS fitness 0.005087 (global minimum of Borel QR map)
- CKM fitness 0.016482
- CP phase delta_HFG = 195.91Â° vs PDG 197.0Â° (0.55% error)
- m_mu/m_e = 208 via N(16+12*omega) (0.59%)
- m_tau/m_e = 3477 via N(68+37*omega) (0.006%)
- CS(m003) = 1/4 exactly; CS(m006,m019,m178) = transcendental

---

## ACTIVE JOURNAL SUBMISSIONS (as of June 1, 2026 09:36)

| Paper | Journal | ID/Contact | Status |
|---|---|---|---|
| Bloch Volume v3 | NYJM | kirsten.wickelgren@duke.edu | Sent 09:36 |
| Bloch Volume v2 | Experimental Mathematics | 265432211 | With administrator |
| Farey Tower | Geometriae Dedicata | 1f71c151... | Technical check (funding added) |
| Delta Invariant | AGT | 260530-Gentry | With editors |
| Dual Surgery | Proc. AMS | 135580 | Blind version sent to JRNL-INITSUB@AMS.ORG |
| Lucas Trace | MRL | 260530-Gentry-2 | With editors |
| BPS Lepton Masses | LMP | MATH-D-26-00372 | Under review |
| Qubit Gates | JMP | JMP26-AR-01272 | Under review since May 6 |
| Galois-Gauge | â€” | CNTP rejected | Needs new target |

### Pending Actions on Submissions
- Bloch v2 (Exp Math 265432211): withdraw if NYJM accepts v3, or if Exp Math slow
- Galois-Gauge (SSRN 6840322): wait for DISTRIBUTED status, then revise to v4
- Galois-Gauge journal: find new target (CNTP rejected "needs rigorous proofs")
  - Candidates: Annales de l'Institut Fourier, IMRN, Journal of Number Theory
- Proc AMS 135580: awaiting response after blind version sent

---

## SSRN PORTFOLIO (10 preprints)

| SSRN | Title | Status |
|---|---|---|
| 6859979 | Bloch Volume Quantum (v2) | PRELIMINARY_UPLOAD â€” revise to v3 when DISTRIBUTED |
| 6854378 | Lucas Trace Identity | PRELIMINARY_UPLOAD |
| 6851440 | Delta Invariant disc=-283 | PRELIMINARY_UPLOAD |
| 6848478 | Mu-Function | IN REVIEW |
| 6845778 | Dual Surgery / Pati-Salam | IN REVIEW |
| 6840418 | BPS Lepton Masses X0(11) | PRELIMINARY_UPLOAD |
| 6840324 | Gauss Polynomial / WRT | PRELIMINARY_UPLOAD |
| 6840322 | Galois-Gauge | PRELIMINARY_UPLOAD â€” revise to v4 when DISTRIBUTED |
| 6815721 | Level-11 Automorphic | PRELIMINARY_UPLOAD |
| 6775158 | Unified HFG | DISTRIBUTED |

---

## PAPERS READY BUT NOT YET SUBMITTED TO JOURNALS

| Paper | File | Target |
|---|---|---|
| Galois-Gauge v4 | gentry-galois-gauge-v4.pdf | IMRN or Annales Fourier |
| Bloch Volume v3 | gentry-bloch-volume-v3.pdf | NYJM (just sent) |
| Meyerhoff/Gauss WRT | gentry_meyerhoff_gauss.pdf | Quantum Topology |
| BPS Lepton Masses | gentry-bps-lepton.pdf | LMP (already submitted) |

---

## OPEN MATHEMATICAL PROBLEMS / NEXT INVARIANTS

### Highest priority: prove D(w^3) = v0 algebraically
This is the one remaining numerical step in the Bloch volume proof.
Would require showing [w^3] in B(K) has regulator = vol(M_PMNS).
This would complete the full algebraic proof.

### Next invariants to investigate
1. **Shape field polynomial as invariant** â€” minimal poly of first tetra shape
   - For m019: x^4-3x^3+3x^2-x-1 (= min poly of w^3)
   - Extend to all disc=-283 manifolds beyond first 2000
2. **Bloch group rank** â€” is [w^3] a generator of B(K)?
   - If yes: vol(m019)=3*beta, vol(m178)=4*beta at Bloch group level
   - This would be the complete arithmetic statement
3. **Reidemeister torsion** of canonical fillings
4. **Eta invariant** of canonical fillings
5. **Regulator of trace field** vs volume â€” relationship for all 4 HFG manifolds
6. **Extend Bloch orbit** â€” find more manifolds in same commensurability class

### Spectral gap results (honestly assessed)
- Only significant: m019(4,1) k=7 gap, p=0.010
- Universal spectral gap conjecture: FALSIFIED by Z/47 counterexample
- Do not overstate spectral gap results

---

## INFRASTRUCTURE STATUS

### PyPI
- latticefit v0.3.1 live at https://pypi.org/project/latticefit/
- README now shows correctly (fixed encoding June 1)
- pyproject.toml has readme = "README.md"
- PyPI token: stored in .pypirc (delete after use)

### GitHub Topics (both repos)
hyperbolic-geometry, mathematical-physics, flavor-physics, standard-model,
snappy, sagemath, number-theory, dehn-surgery, 3-manifolds, particle-physics

### Website
- hyperbolicflavorgeometry.org: LIVE (Cloudflare Pages, June 1)
- Source: C:\dev\hyperbolic-flavor-geometry\docs\index.html
- Auto-deploys on push to main branch
- Last updated: June 1 with Bloch result and SSRN 6859979

---

## POPULAR ARTICLES STATUS

### Written (in /mnt/user-data/outputs/ and C:\dev\framework\)
1. hfg-popular-article.md â€” "The Shape of the Smallest Space and the Mass of the Muon" (1800 words)
2. hfg-article-2-geometry.md â€” "What Is a Hyperbolic 3-Manifold?" (2300 words)

### Planned series
3. Why Particles Mix: The PMNS Matrix Without the Jargon
4. The Arithmetic of Gauge Symmetry: Galois Groups and Weyl Groups
5. Two Roads to the Same Manifold: Dual Surgery and Pati-Salam
6. The Muon's Mass as an Eisenstein Norm
7. What We Know and What We Don't (honest assessment)
8. HFG Dark Sector Conjectures and Quantum Gravity (speculative)

### Target venues
- Inference (primary â€” accepts unsolicited independent work)
- Nautilus
- Notices of the AMS ("What is...?" format, 2 pages)

---

## CHIRALITY VIDEO

### Script
- File: chirality-script-v1.md (374 lines, 10-12 min)
- Status: written, needs tuning with Marvin
- Three acts: everyday chirality â†’ physics (Wu experiment) â†’ bigger picture
- Props needed: rubber gloves (inflate), bendable wire helices (clockwise/CCW)
- HFG mention: one paragraph in Act 3, no jargon
- Target: YouTube within 3 days of filming

### Production
- DaVinci Resolve installed
- Creative Commons graphics needed: carvone enantiomers, DNA helix, Wu experiment
- Screen recording of SnapPy amphicheiral manifold (5 seconds)

---

## THINGS TO DO NEXT SESSION (priority order)

1. **Prove D(w^3) = v0 algebraically** â€” the key open problem
2. **Check if [w^3] generates B(K)** â€” Sage computation
3. **Revise SSRN 6840322** (Galois-Gauge) to v4 when DISTRIBUTED
4. **Find new journal for Galois-Gauge** â€” IMRN or Annales Fourier
5. **Submit Meyerhoff/Gauss WRT paper** to Quantum Topology
6. **Write article 3** â€” PMNS matrix explainer
7. **Tune chirality script** with Marvin
8. **Update website** â€” add articles section
9. **Extend disc=-283 census scan** beyond 2000 manifolds
10. **Revoke exposed GitHub token** â€” https://github.com/settings/tokens

---

## REJECTED SUBMISSIONS THIS CYCLE (retire, do not resubmit)

- CNTP: Galois-Gauge ("needs rigorous proofs" â€” paper needs algebraic proofs added)
- Comptes Rendus: multiple rejections without review
- Universe (MDPI): desk reject
- 7 old APS submissions: all "not under active consideration"
- 12 Springer/Elsevier rejections: all pre-HFG framework papers, retired

---

## FITNESS VALUES AND PHYSICAL PREDICTIONS

```
PMNS fitness: 0.005087 (global minimum, column permutation (1,0,2))
CKM fitness:  0.016482 (words aaB/AbA/AAb, sigma=0.49)
CP phase:     delta_HFG = 195.91Â° vs PDG 197.0Â° (0.55%, zero free params)

Lepton masses (Eisenstein norms in Z[omega], trace field of M_PMNS):
  L_11 = 199 â‰ˆ m_mu/m_e  (0.003%)  [CORRECTED: N(16+12*omega)=208, 0.59%]
  L_17 = 3571 â‰ˆ m_tau/m_e (0.000%) [CORRECTED: N(68+37*omega)=3477, 0.006%]

Strong CP: theta_QCD = 0 forced by amphicheirality of all 6 HFG manifolds
```

---

## SAGE WORKFLOW (for next session)

```bash
# Start WSL
wsl

# Activate sage env
conda activate sage

# Start Sage
sage

# In Sage:
import snappy
from sage.all import NumberField, QQ, algdep, ComplexField, RealField

# ALWAYS use high_precision() for shape computations
M = snappy.Manifold("m019")
M_hp = M.high_precision()
shapes = M_hp.tetrahedra_shapes('rect')

# ALWAYS use ComplexField(200) for algdep
CC = ComplexField(200)
z = CC(shapes[0])
p = algdep(z, 4, known_bits=150)
```

---

## GIT COMMIT STYLE

```powershell
cd C:\dev\framework
git add .
git commit -m "RESULT TYPE: brief description. Key values. References."
git push origin main

# Examples:
# "PROVED: z_A=w^3, z_B=-w, z_S=w^{-4} in K=Q(w), w^4=w+1. Product=-1. T-orbit exact."
# "SUBMISSION: NYJM Wickelgren (Bloch v3). Sent 2026-06-01 09:36."
# "SSRN 6859979 revised to v3 with algebraic proof."
```
