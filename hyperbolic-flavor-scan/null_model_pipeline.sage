# null_model_pipeline.sage
# DS-designed statistical pipeline for arithmetic norm matching
# Replaces fragile pattern-matching with proper null model comparison
# Run: sage null_model_pipeline.sage

from sage.all import *
import random
import numpy as np

# ============================================================
# 0. DATA: fermion mass ratios (m / m_e)
# ============================================================
MASSES = {
    "mu":  206.77,
    "tau": 3477.22,
    "u":   4.24,
    "d":   9.14,
    "s":   182.78,
    "c":   2495.11,
    "b":   8180.04,
    "t":   338082.19,
}

# ============================================================
# 1. SAFE NORM DEFINITIONS (no symbolic expressions)
# ============================================================
def norm_Zsqrt_d(a, b, d):
    """Real quadratic norm: |a^2 - d*b^2|"""
    return abs(int(a)*int(a) - int(d)*int(b)*int(b))

def norm_eisenstein(a, b):
    """Eisenstein norm: a^2 - ab + b^2"""
    a,b = int(a),int(b)
    return a*a - a*b + b*b

def norm_gaussian(a, b):
    """Gaussian norm: a^2 + b^2"""
    return int(a)*int(a) + int(b)*int(b)

# ============================================================
# 2. WITNESS SEARCH (bounded, deterministic)
# ============================================================
def find_best_witness(m, norm_func, R=100):
    """Brute-force bounded lattice search."""
    best_err = 1e100
    best = (0, 0)
    best_n = 0
    m = float(m)
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            n = norm_func(a, b)
            if n <= 0:
                continue
            err = abs(n - m) / m
            if err < best_err:
                best_err = err
                best = (a, b)
                best_n = n
    return best_n, best, best_err

# ============================================================
# 3. STRUCTURED MODELS
# ============================================================
STRUCTURED_MODELS = {
    "Z[omega] (Eisenstein)": lambda a,b: norm_eisenstein(a,b),
    "Z[sqrt(17)] (CKM)":     lambda a,b: norm_Zsqrt_d(a,b,17),
    "Z[sqrt(5)] (phi)":      lambda a,b: norm_Zsqrt_d(a,b,5),
    "Z[i] (Gaussian)":       lambda a,b: norm_gaussian(a,b),
}

# ============================================================
# 4. NULL MODELS
# ============================================================
def random_quadratic_factory(d_val):
    """Norm form |a^2 - d*b^2| for random d (not 17 or 5 or 3)."""
    def model(a, b):
        return norm_Zsqrt_d(a, b, d_val)
    return model

NULL_MODELS = {
    "Z[sqrt(2)]":   random_quadratic_factory(2),
    "Z[sqrt(3)]":   random_quadratic_factory(3),
    "Z[sqrt(7)]":   random_quadratic_factory(7),
    "Z[sqrt(11)]":  random_quadratic_factory(11),
    "Z[sqrt(13)]":  random_quadratic_factory(13),
    "Z[sqrt(19)]":  random_quadratic_factory(19),
    "Z[sqrt(23)]":  random_quadratic_factory(23),
    "Z[sqrt(29)]":  random_quadratic_factory(29),
    "Z[sqrt(31)]":  random_quadratic_factory(31),
    "Z[sqrt(37)]":  random_quadratic_factory(37),
}

# ============================================================
# 5. EXPERIMENT RUNNER
# ============================================================
def run_model(name, model_func, masses, R=100, verbose=False):
    results = {}
    if verbose:
        print(f"\n{name}:")
    for f, m in masses.items():
        n, (a,b), err = find_best_witness(m, model_func, R=R)
        results[f] = {"mass": float(m), "norm": n,
                      "witness": (a,b), "error": float(err)}
        if verbose:
            print(f"  {f:>3}: m={m:10.2f}  N={n:10}  err={err:.4%}")
    return results

def mean_log_error(results):
    errors = [v["error"] for v in results.values()]
    return float(np.mean(np.log10(errors)))

def median_error(results):
    errors = [v["error"] for v in results.values()]
    return float(np.median(errors))

def fraction_within(results, tol):
    errors = [v["error"] for v in results.values()]
    return sum(1 for e in errors if e < tol) / len(errors)

# ============================================================
# 6. MAIN PIPELINE
# ============================================================
print("="*70)
print("NULL MODEL PIPELINE: Arithmetic Norm Matching for Fermion Masses")
print("="*70)

# ── Run structured models ────────────────────────────────────────────────────
print("\nSTRUCTURED MODELS (detailed):")
structured_stats = {}
for name, model in STRUCTURED_MODELS.items():
    res = run_model(name, model, MASSES, R=100, verbose=True)
    mle = mean_log_error(res)
    med = median_error(res)
    f2 = fraction_within(res, 0.02)
    f5 = fraction_within(res, 0.05)
    structured_stats[name] = (mle, med, f2, f5, res)

# ── Run null models ──────────────────────────────────────────────────────────
print("\nNULL MODELS (random quadratic forms):")
null_stats = {}
for name, model in NULL_MODELS.items():
    res = run_model(name, model, MASSES, R=100, verbose=False)
    mle = mean_log_error(res)
    med = median_error(res)
    f2 = fraction_within(res, 0.02)
    f5 = fraction_within(res, 0.05)
    null_stats[name] = (mle, med, f2, f5, res)
    print(f"  {name:<18}: mean_log_err={mle:.3f}  median={med:.4f}  "
          f"frac<2%={f2:.2f}  frac<5%={f5:.2f}")

# ── Shuffled mass null ───────────────────────────────────────────────────────
print("\nSHUFFLED MASS NULL (permutation test, N=1000):")
observed_f2 = {name: stats[2] for name, stats in structured_stats.items()}
observed_mle = {name: stats[0] for name, stats in structured_stats.items()}

N_SHUFFLE = 1000
shuffle_f2 = {name: [] for name in STRUCTURED_MODELS}
shuffle_mle = {name: [] for name in STRUCTURED_MODELS}

masses_list = list(MASSES.values())
keys_list = list(MASSES.keys())

for trial in range(N_SHUFFLE):
    random.seed(trial)
    shuffled_vals = masses_list[:]
    random.shuffle(shuffled_vals)
    shuffled = dict(zip(keys_list, shuffled_vals))
    for name, model in STRUCTURED_MODELS.items():
        res = run_model(name, model, shuffled, R=60, verbose=False)
        shuffle_f2[name].append(fraction_within(res, 0.02))
        shuffle_mle[name].append(mean_log_error(res))

print(f"\n{'Model':<30} {'Obs f<2%':>9} {'Null mean':>9} "
      f"{'p-value':>8} {'Obs MLE':>9} {'Null MLE':>9} {'p-value':>8}")
print("-"*80)
for name in STRUCTURED_MODELS:
    obs_f = observed_f2[name]
    null_mean_f = float(np.mean(shuffle_f2[name]))
    p_f = sum(1 for x in shuffle_f2[name] if x >= obs_f) / N_SHUFFLE
    obs_m = observed_mle[name]
    null_mean_m = float(np.mean(shuffle_mle[name]))
    p_m = sum(1 for x in shuffle_mle[name] if x <= obs_m) / N_SHUFFLE
    print(f"{name:<30} {obs_f:>9.2f} {null_mean_f:>9.2f} "
          f"{p_f:>8.4f} {obs_m:>9.3f} {null_mean_m:>9.3f} {p_m:>8.4f}")

# ── Summary table ────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("SUMMARY: STRUCTURED vs NULL MODELS")
print("="*70)
null_mle_vals = [s[0] for s in null_stats.values()]
null_f2_vals = [s[2] for s in null_stats.values()]
null_mean_mle = float(np.mean(null_mle_vals))
null_mean_f2 = float(np.mean(null_f2_vals))

print(f"\nNull model average: mean_log_err={null_mean_mle:.3f}, frac<2%={null_mean_f2:.3f}")
print()
print(f"{'Model':<30} {'MLE':>8} {'vs null':>8} {'f<2%':>7} {'vs null':>8}")
print("-"*65)
for name, (mle, med, f2, f5, res) in structured_stats.items():
    delta_mle = mle - null_mean_mle
    delta_f2 = f2 - null_mean_f2
    better = "BETTER" if delta_mle < 0 else "worse"
    print(f"{name:<30} {mle:>8.3f} {delta_mle:>+8.3f} {f2:>7.2f} "
          f"{delta_f2:>+8.2f}  {better}")

print()
print("="*70)
print("LATEX TABLE: Best witnesses per fermion")
print("="*70)
print()
print(r"\begin{tabular}{llrrrr}")
print(r"\toprule")
print(r"Fermion & Ring & $m/m_e$ & Norm $N$ & Witness & Error \\")
print(r"\midrule")

# Leptons use Eisenstein, quarks use Z[sqrt(17)]
best_lepton = structured_stats["Z[omega] (Eisenstein)"][4]
best_quark  = structured_stats["Z[sqrt(17)] (CKM)"][4]

for f in ["mu", "tau"]:
    v = best_lepton[f]
    a,b = v['witness']
    print(f"{f} & $\\mathbb{{Z}}[\\omega]$ & {v['mass']:.2f} & "
          f"{v['norm']} & $({a},{b})$ & {v['error']*100:.4f}\\% \\\\")
for f in ["u","d","s","c","b","t"]:
    v = best_quark[f]
    a,b = v['witness']
    print(f"{f} & $\\mathbb{{Z}}[\\sqrt{{17}}]$ & {v['mass']:.2f} & "
          f"{v['norm']} & $({a},{b})$ & {v['error']*100:.4f}\\% \\\\")
print(r"\bottomrule")
print(r"\end{tabular}")
