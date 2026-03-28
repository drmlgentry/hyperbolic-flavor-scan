"""
stronger_null_tests.py
======================
Runs three null tests that are missing from the current papers:

1. Haar-random unitary null (CKM and PMNS)
   - Generates random unitary matrices from the Haar measure
   - Tests whether the geometric fit is better than random unitaries
   - This is the physically motivated null for a mixing matrix claim

2. Log-normal null with matched moments (mass lattice)
   - Draws masses from a log-normal with same mean/variance as SM data
   - Stronger than log-uniform — tests against a realistic prior

3. Bootstrap on PDG error bars (mass lattice)
   - Resamples masses within their experimental uncertainties
   - Tests stability of the lattice fit to measurement errors

4. Trials correction (base/denominator search)
   - Documents every base and denominator tested
   - Applies Bonferroni correction to reported p-values

Usage:
    cd C:\\dev\\hyperbolic-flavor-scan
    python analysis\\stronger_null_tests.py

Output:
    analysis\\stronger_null_results.txt
    analysis\\stronger_null_results.csv
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import unitary_group
from scipy.linalg import qr
import warnings
warnings.filterwarnings('ignore')

rng = np.random.default_rng(42)
PHI = (1 + 5**0.5) / 2

OUT_DIR = Path(__file__).parent
RESULTS = []

def log(*args):
    msg = " ".join(str(a) for a in args)
    print(msg)
    return msg


# ── PDG 2024 values ───────────────────────────────────────────────
# CKM matrix (standard parametrisation, PDG 2024)
CKM_PDG = np.array([
    [0.97435, 0.22501, 0.00369],
    [0.22487, 0.97350, 0.04182],
    [0.00857, 0.04110, 0.99912],
], dtype=float)  # magnitudes |V_ij|

# PMNS matrix (magnitudes, PDG 2024 best fit)
PMNS_PDG = np.array([
    [0.8249, 0.5479, 0.1521],
    [0.4685, 0.5756, 0.6714],
    [0.3152, 0.6070, 0.7259],
], dtype=float)

# SM fermion masses (GeV, PDG 2024) with uncertainties
SM_MASSES = np.array([
    5.11e-4, 0.10566, 1.77686,   # e, mu, tau
    0.00216, 1.275,   172.76,    # u, c, t
    0.00467, 0.0934,  4.18,      # d, s, b
])
# Fractional uncertainties (approximate, from PDG)
SM_MASS_ERR = np.array([
    0.0,    0.0,    0.00006,   # e (exact), mu (exact), tau
    0.25,   0.04,   0.0009,    # u (large), c, t
    0.15,   0.06,   0.008,     # d (large), s, b
])

# Predicted lattice labels q (from the paper)
SM_Q = np.array([0, 44, 68, 12, 65, 106, 18, 43, 75])


# ── Helper: matrix fitness ────────────────────────────────────────
def matrix_fitness(M_pred, M_obs):
    """RMS difference between predicted and observed matrix magnitudes."""
    return float(np.sqrt(np.mean((np.abs(M_pred) - np.abs(M_obs))**2)))


# ── Helper: lattice RMS ───────────────────────────────────────────
def lattice_rms(vals, base, denom, anchor):
    """RMS log-residual of vals from the lattice anchor*base^(k/denom)."""
    vals = np.asarray(vals, dtype=float)
    vals = vals[vals > 0]
    ks   = np.round(denom * np.log(vals/anchor) / np.log(base))
    pred = anchor * base**(ks/denom)
    res  = np.abs(np.log(vals/pred) / np.log(base))
    return float(np.sqrt(np.mean(res**2)))


# ── 1. HAAR-RANDOM UNITARY NULL (CKM) ────────────────────────────
log("\n" + "="*60)
log("TEST 1: Haar-random unitary null — CKM matrix")
log("="*60)
log("Question: Is the geometric fit to CKM better than random")
log("          unitary matrices drawn from the Haar measure?")
log()

# Geometric CKM prediction (from the paper, best word triple)
# Using the reported fitness value from the CKM paper
CKM_GEOMETRIC_FITNESS = 0.01729  # from word_triple_results_corrected.csv

# Generate Haar-random 3x3 unitary matrices
n_trials = 50000
haar_fitness = []
for _ in range(n_trials):
    # Draw from Haar measure via QR decomposition of complex Gaussian
    Z = (rng.standard_normal((3,3)) + 1j*rng.standard_normal((3,3))) / np.sqrt(2)
    Q, R = qr(Z)
    # Adjust phase to get uniform Haar measure
    D = np.diag(R.diagonal() / np.abs(R.diagonal()))
    U = Q @ D
    haar_fitness.append(matrix_fitness(np.abs(U), CKM_PDG))

haar_fitness = np.array(haar_fitness)
p_haar_ckm = float(np.mean(haar_fitness <= CKM_GEOMETRIC_FITNESS))

log(f"Geometric fitness (CKM): {CKM_GEOMETRIC_FITNESS:.5f}")
log(f"Haar null: mean={haar_fitness.mean():.5f}, "
    f"std={haar_fitness.std():.5f}, "
    f"min={haar_fitness.min():.5f}")
log(f"p-value (Haar null, CKM): {p_haar_ckm:.4f}")
if p_haar_ckm < 0.001:
    log("  >>> HIGHLY SIGNIFICANT vs Haar measure")
elif p_haar_ckm < 0.05:
    log("  >>> SIGNIFICANT vs Haar measure")
else:
    log("  >>> NOT significant vs Haar measure")

RESULTS.append({
    "test": "Haar-random unitary (CKM)",
    "observed_stat": CKM_GEOMETRIC_FITNESS,
    "null_mean": float(haar_fitness.mean()),
    "null_std": float(haar_fitness.std()),
    "p_value": p_haar_ckm,
    "n_trials": n_trials,
    "significant_0.05": p_haar_ckm < 0.05,
})

# ── 1b. HAAR-RANDOM UNITARY NULL (PMNS) ──────────────────────────
log()
log("-"*40)
log("TEST 1b: Haar-random unitary null — PMNS matrix")
log("-"*40)

PMNS_GEOMETRIC_FITNESS = 0.01897  # from the PMNS paper

pmns_haar_fitness = []
for _ in range(n_trials):
    Z = (rng.standard_normal((3,3)) + 1j*rng.standard_normal((3,3))) / np.sqrt(2)
    Q, R = qr(Z)
    D = np.diag(R.diagonal() / np.abs(R.diagonal()))
    U = Q @ D
    pmns_haar_fitness.append(matrix_fitness(np.abs(U), PMNS_PDG))

pmns_haar_fitness = np.array(pmns_haar_fitness)
p_haar_pmns = float(np.mean(pmns_haar_fitness <= PMNS_GEOMETRIC_FITNESS))

log(f"Geometric fitness (PMNS): {PMNS_GEOMETRIC_FITNESS:.5f}")
log(f"Haar null: mean={pmns_haar_fitness.mean():.5f}, "
    f"std={pmns_haar_fitness.std():.5f}")
log(f"p-value (Haar null, PMNS): {p_haar_pmns:.4f}")
if p_haar_pmns < 0.001:
    log("  >>> HIGHLY SIGNIFICANT vs Haar measure")
elif p_haar_pmns < 0.05:
    log("  >>> SIGNIFICANT vs Haar measure")
else:
    log("  >>> NOT significant vs Haar measure")

RESULTS.append({
    "test": "Haar-random unitary (PMNS)",
    "observed_stat": PMNS_GEOMETRIC_FITNESS,
    "null_mean": float(pmns_haar_fitness.mean()),
    "null_std": float(pmns_haar_fitness.std()),
    "p_value": p_haar_pmns,
    "n_trials": n_trials,
    "significant_0.05": p_haar_pmns < 0.05,
})

# ── 2. LOG-NORMAL NULL (mass lattice) ─────────────────────────────
log()
log("="*60)
log("TEST 2: Log-normal null with matched moments — mass lattice")
log("="*60)
log("Question: Is the lattice fit better than masses drawn from")
log("          a log-normal with the same mean and variance?")
log("          (Stronger than log-uniform — physically motivated)")
log()

masses = SM_MASSES
anchor = masses[0]  # electron mass as anchor
base   = PHI
denom  = 4

obs_rms = lattice_rms(masses, base, denom, anchor)

# Fit log-normal to the observed masses
log_masses = np.log(masses)
mu_ln  = log_masses.mean()
sig_ln = log_masses.std()

log(f"Log-normal fit: mu={mu_ln:.4f}, sigma={sig_ln:.4f}")
log(f"Observed RMS (phi lattice): {obs_rms:.5f}")

lognorm_rms = []
for _ in range(n_trials):
    # Draw 9 masses from log-normal with matched moments
    r = np.exp(rng.normal(mu_ln, sig_ln, len(masses)))
    r = np.sort(r)  # sort to match ordering
    lognorm_rms.append(lattice_rms(r, base, denom, r[0]))

lognorm_rms = np.array(lognorm_rms)
p_lognorm = float(np.mean(lognorm_rms <= obs_rms))

log(f"Log-normal null: mean={lognorm_rms.mean():.5f}, "
    f"std={lognorm_rms.std():.5f}")
log(f"p-value (log-normal null): {p_lognorm:.4f}")
if p_lognorm < 0.001:
    log("  >>> HIGHLY SIGNIFICANT vs log-normal")
elif p_lognorm < 0.05:
    log("  >>> SIGNIFICANT vs log-normal")
elif p_lognorm < 0.10:
    log("  >>> MARGINAL vs log-normal")
else:
    log("  >>> NOT significant vs log-normal")

RESULTS.append({
    "test": "Log-normal null (mass lattice)",
    "observed_stat": obs_rms,
    "null_mean": float(lognorm_rms.mean()),
    "null_std": float(lognorm_rms.std()),
    "p_value": p_lognorm,
    "n_trials": n_trials,
    "significant_0.05": p_lognorm < 0.05,
})

# ── 3. PDG BOOTSTRAP NULL (mass lattice) ─────────────────────────
log()
log("="*60)
log("TEST 3: Bootstrap on PDG error bars — mass lattice")
log("="*60)
log("Question: Is the lattice fit stable under measurement")
log("          uncertainties? Does the p-value hold if masses")
log("          are perturbed within their PDG errors?")
log()

bootstrap_rms  = []
bootstrap_p    = []
n_boot = 10000

for _ in range(n_boot):
    # Perturb each mass within its fractional uncertainty
    perturbed = masses * np.exp(
        rng.normal(0, SM_MASS_ERR, len(masses))
    )
    perturbed = perturbed[perturbed > 0]
    if len(perturbed) < 5:
        continue
    rms_b = lattice_rms(perturbed, base, denom, perturbed[0])
    bootstrap_rms.append(rms_b)

bootstrap_rms = np.array(bootstrap_rms)
log(f"Bootstrap RMS: mean={bootstrap_rms.mean():.5f} "
    f"± {bootstrap_rms.std():.5f}")
log(f"95% CI: [{np.percentile(bootstrap_rms,2.5):.5f}, "
    f"{np.percentile(bootstrap_rms,97.5):.5f}]")
log(f"Observed RMS: {obs_rms:.5f}")
log(f"Fraction of bootstrap samples with RMS < observed: "
    f"{np.mean(bootstrap_rms < obs_rms):.3f}")
log(f"Result: lattice fit is {'STABLE' if bootstrap_rms.std() < 0.01 else 'UNSTABLE'} "
    f"under PDG uncertainties")

RESULTS.append({
    "test": "PDG bootstrap (mass lattice)",
    "observed_stat": obs_rms,
    "null_mean": float(bootstrap_rms.mean()),
    "null_std": float(bootstrap_rms.std()),
    "p_value": float(np.mean(bootstrap_rms < obs_rms)),
    "n_trials": n_boot,
    "significant_0.05": bootstrap_rms.std() < 0.005,
})

# ── 4. TRIALS CORRECTION ─────────────────────────────────────────
log()
log("="*60)
log("TEST 4: Trials correction for base/denominator search")
log("="*60)
log("Documenting all bases and denominators tested,")
log("applying Bonferroni correction to reported p-value.")
log()

BASES_TESTED  = [PHI, 2**0.5, 2.0, np.e, np.pi, 10.0, PHI**2]
DENOMS_TESTED = [2, 3, 4, 6, 8, 12]
N_TRIALS_SEARCH = len(BASES_TESTED) * len(DENOMS_TESTED)

log(f"Bases tested: {len(BASES_TESTED)} "
    f"(phi, sqrt2, 2, e, pi, 10, phi^2)")
log(f"Denominators tested: {len(DENOMS_TESTED)} (2,3,4,6,8,12)")
log(f"Total combinations: {N_TRIALS_SEARCH}")

# Run all combinations on SM masses
all_results = []
for b in BASES_TESTED:
    for d in DENOMS_TESTED:
        rms = lattice_rms(masses, b, d, anchor)
        all_results.append((b, d, rms))

all_results.sort(key=lambda x: x[2])
log(f"\nTop 5 combinations:")
for b, d, rms in all_results[:5]:
    bname = {PHI:"phi", 2**0.5:"sqrt2", 2.0:"2",
             np.e:"e", np.pi:"pi", 10.0:"10", PHI**2:"phi^2"}.get(b, f"{b:.4f}")
    log(f"  base={bname}, d={d}, RMS={rms:.5f}")

# Best observed RMS
best_rms = all_results[0][2]
best_b   = all_results[0][0]
best_d   = all_results[0][1]

# Uncorrected p-value (log-uniform null for best combo)
lmin, lmax = np.log(masses.min()), np.log(masses.max())
uncorr_null = []
for _ in range(n_trials):
    r = np.exp(rng.uniform(lmin, lmax, len(masses)))
    uncorr_null.append(lattice_rms(r, best_b, best_d, r[0]))
p_uncorrected = float(np.mean(np.array(uncorr_null) <= best_rms))

# Bonferroni corrected
p_bonferroni = min(p_uncorrected * N_TRIALS_SEARCH, 1.0)

# Holm-Bonferroni (less conservative)
# Compute p for each combination
combo_p = []
for b, d, rms in all_results:
    null_i = []
    for _ in range(5000):
        r = np.exp(rng.uniform(lmin, lmax, len(masses)))
        null_i.append(lattice_rms(r, b, d, r[0]))
    combo_p.append(float(np.mean(np.array(null_i) <= rms)))

# Holm step-down
combo_p_sorted = sorted(enumerate(combo_p), key=lambda x: x[1])
holm_p = []
n = len(combo_p)
for rank, (orig_idx, p) in enumerate(combo_p_sorted):
    holm_p.append(min(p * (n - rank), 1.0))
p_holm_best = holm_p[0]

log(f"\nReported p-value (log-uniform, best combo): {p_uncorrected:.4f}")
log(f"Bonferroni corrected ({N_TRIALS_SEARCH} trials): {p_bonferroni:.4f}")
log(f"Holm-Bonferroni corrected: {p_holm_best:.4f}")

if p_bonferroni < 0.05:
    log("  >>> SIGNIFICANT after Bonferroni correction")
elif p_holm_best < 0.05:
    log("  >>> SIGNIFICANT after Holm correction (but not Bonferroni)")
else:
    log("  >>> NOT significant after trials correction")
    log("  NOTE: This is the most important finding.")
    log("  The reported p-values in the papers do not account for")
    log("  the search over bases and denominators.")

RESULTS.append({
    "test": "Trials correction (Bonferroni)",
    "observed_stat": best_rms,
    "null_mean": float(np.mean(uncorr_null)),
    "null_std": float(np.std(uncorr_null)),
    "p_value": p_bonferroni,
    "n_trials": N_TRIALS_SEARCH,
    "significant_0.05": p_bonferroni < 0.05,
})
RESULTS.append({
    "test": "Trials correction (Holm)",
    "observed_stat": best_rms,
    "null_mean": float(np.mean(uncorr_null)),
    "null_std": float(np.std(uncorr_null)),
    "p_value": p_holm_best,
    "n_trials": N_TRIALS_SEARCH,
    "significant_0.05": p_holm_best < 0.05,
})

# ── SUMMARY ───────────────────────────────────────────────────────
log()
log("="*60)
log("SUMMARY OF ALL NULL TESTS")
log("="*60)
log(f"{'Test':<35} {'p-value':>9} {'Sig?':>6}")
log("-"*55)
for r in RESULTS:
    sig = "YES" if r['significant_0.05'] else "no"
    log(f"  {r['test']:<33} {r['p_value']:>9.4f} {sig:>6}")

# Save results
df = pd.DataFrame(RESULTS)
df.to_csv(OUT_DIR / "stronger_null_results.csv", index=False)

# Save text report
report_lines = []
report_lines.append("STRONGER NULL TEST RESULTS")
report_lines.append("Hyperbolic Flavor Geometry framework")
report_lines.append(f"Date: March 2026 · n_trials={n_trials}")
report_lines.append("="*60)
for r in RESULTS:
    report_lines.append(
        f"{r['test']}: p={r['p_value']:.4f} "
        f"({'SIGNIFICANT' if r['significant_0.05'] else 'not significant'})"
    )
report_lines.append("")
report_lines.append("INTERPRETATION:")
report_lines.append(
    "The Haar-random unitary null tests whether the mixing matrix")
report_lines.append(
    "fitness is better than a random unitary. If p < 0.05 here,")
report_lines.append(
    "the result is genuinely non-trivial. If not, the geometric")
report_lines.append(
    "construction does not outperform random unitaries.")
report_lines.append("")
report_lines.append(
    "The trials correction is critical for the mass lattice claim.")
report_lines.append(
    "If p_Bonferroni > 0.05, the reported p-values in the papers")
report_lines.append(
    "are not corrected for the search and should be revised.")

with open(OUT_DIR / "stronger_null_results.txt", "w") as f:
    f.write("\n".join(report_lines))

log(f"\nResults saved to:")
log(f"  {OUT_DIR}/stronger_null_results.csv")
log(f"  {OUT_DIR}/stronger_null_results.txt")
