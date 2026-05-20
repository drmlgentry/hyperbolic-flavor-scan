# Reproducibility Guide — Farey Tower Paper
**Paper:** "A Quadratic Cover Prime Formula for a Farey Tower of Arithmetic Hyperbolic 3-Manifolds"  
**Author:** Marvin L. Gentry  
**File:** `paper/gentry-farey-tower.tex`  
**GitHub:** https://github.com/drmlgentry/hyperbolic-flavor-scan

---

## Environment

All scripts require SnapPy and SageMath in a conda environment.

```bash
# Activate environment
conda activate sage

# Verify SnapPy
conda run -n sage python -c "import snappy; print(snappy.__version__)"

# Working directory (WSL path)
cd /mnt/c/dev/hyperbolic-flavor-scan
```

SnapPy version used: 3.3.2  
Python: 3.x in conda sage environment  
Platform: WSL2 on Windows (Ubuntu 24)

---

## Scripts and What They Prove

### 1. `cover_prime_verify.py` — Main verification (canonical)

**Verifies Theorem 6.1** (the cover homology formula):
```
H1(M_n^(5)) = (Z/p_n)^2,  p_n = 5*n*(n+1) + 1
```
for n = 1..6 (default). Checks:
- Each M_n = m003(-n, 2n+1) is hyperbolic with H1 = Z/5
- No algebraic covers of degree 2, 3, or 4 exist
- Exactly one degree-5 algebraic cover exists
- Cover volume = 5 * vol(M_n) exactly (genuine covering relation)
- Cover H1 = (Z/p_n)^2 matching the formula
- Lucas-primality: only n=1 gives p_1=11=L_5 (Lucas prime)

```bash
conda run -n sage python cover_prime_verify.py 2>/dev/null
```

**Expected output (last lines):**
```
ALL 6 CASES PASS.
Formula p_n=5n(n+1)+1 verified for n=1..6.
Lucas-primality: ONLY n=1 (p_1=11=L_5).
```

**Runtime:** ~10 seconds  
**Output file:** `data/cover_prime_verification.csv`

---

### 2. `linking_form.py` — Surgery formula proof

**Verifies and proves Theorem 4.1** (surgery formula):
```
|H1(m003(p,q))| = 5 * |2p + q|
```

Derives the formula algebraically from:
- Fundamental group: generators a,b; relator abA^2babbb
- Meridian word: ABABB → abelianizes to -2α - 3β
- Longitude word: ABAbab → abelianizes to -α + β
- Surgery determinant: |det[(0,5),(-2p-q,-3p+q)]| = 5|2p+q|

Also verifies against 14 independent fillings (all match).

```bash
conda run -n sage python linking_form.py 2>/dev/null
```

**Expected output (last lines):**
```
ALL SLOPES MATCH. Formula |H1| = 5|2p+q| is exact.
...
Corollary: H1(m003(-n,2n+1)) = Z/5 for ALL n >= 1.
```

**Runtime:** ~5 seconds

---

### 3. `farey_cover_test.py` — Extended cover prime test

Tests the cover prime formula across multiple slopes of m003,
including the Farey chain and neighboring slopes. Also tests
the Farey intersection conjecture (det[slope1, slope2] = cover prime).

```bash
conda run -n sage python farey_cover_test.py 2>/dev/null
```

**Output file:** `data/farey_cover_results.csv`  
**Runtime:** ~5 minutes (cover enumeration to degree 6)

---

### 4. `h1_slope_formula.py` — Surgery formula landscape

Maps H1(m003(p,q)) across all coprime slopes |p|,|q| <= 7.
Identifies the six H1=Z/5 slopes and their volumes.
Verifies the three distinct isometry classes (Farey chain structure).

```bash
conda run -n sage python h1_slope_formula.py 2>/dev/null
```

**Output file:** `data/h1_slope_results.csv`  
**Runtime:** ~30 seconds

---

## Extended Verification (n=1..10)

The formula was additionally verified for n=7..10 using:

```bash
python3 -c "
import snappy
from sympy import factorint

def full_cover_h1(n):
    M = snappy.Manifold('m003')
    M.dehn_fill((-n, 2*n+1))
    base_order = M.homology().order()
    covers = M.covers(5)
    return [(str(c.homology()), c.homology().order()) for c in covers]

for n in range(1, 11):
    pn = 5*n*(n+1)+1
    info = full_cover_h1(n)
    print(f'n={n:2d} p_n={pn:5d}={factorint(pn)}  covers={info}')
" 2>/dev/null
```

**Result:** H1 = (Z/p_n)^2 for all n=1..10, including composite p_n.

---

## Key Numerical Results (for paper cross-check)

| Claim | Value | Script |
|---|---|---|
| Cusp shape tau | 0.5 + 0.866025...i, residual 9.55e-16 | `linking_form.py` |
| Cusp translations | mu=1+sqrt(3)*i, lam=2 | `linking_form.py` |
| Cusp area | 2*sqrt(3) = 3.4641..., residual 3.11e-15 | `linking_form.py` |
| Surgery formula | 5*|2p+q|, 14/14 slopes match | `linking_form.py` |
| Cover n=1: p_1 | 11 = L_5 (Lucas prime) | `cover_prime_verify.py` |
| Cover n=6: p_6 | 211 (prime) | `cover_prime_verify.py` |
| Cover n=8: p_8 | 361 = 19^2, H1=(Z/361)^2 | extended inline script |
| Cover n=9: p_9 | 451 = 11x41, H1=(Z/451)^2 | extended inline script |
| vol(M_1) | 0.98136883 | any script |
| vol(m003 cusp) | 2.02988321 | `h1_slope_formula.py` |
| NZ constant C | ~14.16 (gap*L^2 converging) | inline computation |

---

## Data Files

```
data/
  cover_prime_verification.csv   -- main theorem evidence (n=1..6)
  farey_cover_results.csv        -- Farey/cover conjecture test
  h1_slope_results.csv           -- surgery formula landscape
  FAREY_TOWER_THEOREM.md         -- theorem documentation
```

---

## Compiling the Paper

```bash
# In any LaTeX environment with standard packages
pdflatex gentry-farey-tower.tex
pdflatex gentry-farey-tower.tex  # second pass for cross-references
```

Required packages: amsmath, amssymb, amsthm, authblk, geometry,
hyperref, booktabs, enumerate, xcolor. All standard in TeX Live 2020+.

---

## Contact

Marvin L. Gentry  
drmlgentry@protonmail.com  
ORCID: 0009-0006-4550-2663
