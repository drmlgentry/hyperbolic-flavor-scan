# null_model_pipeline_v2.sage
# Fixed: force all arithmetic to Python int/float, not Sage types
# Extended radius for heavy fermions
# Run: sage null_model_pipeline_v2.sage

from sage.all import *
import random, numpy as np

# ── Force Python types throughout ────────────────────────────────────────────
def to_float(x):
    try: return float(x)
    except: return float(str(x))

def to_int(x):
    try: return int(x)
    except: return int(str(x))

# ── Mass ratios ───────────────────────────────────────────────────────────────
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

# ── Norm functions (Python int arithmetic only) ───────────────────────────────
def norm_sqrt_d(a, b, d):
    a,b,d = to_int(a), to_int(b), to_int(d)
    return abs(a*a - d*b*b)

def norm_eisenstein(a, b):
    a,b = to_int(a), to_int(b)
    n = a*a - a*b + b*b
    return n if n > 0 else 0

def norm_gaussian(a, b):
    a,b = to_int(a), to_int(b)
    return a*a + b*b

# ── Witness search ────────────────────────────────────────────────────────────
def find_best(m, norm_func, R=150):
    m = float(m)
    best_err = 1e100
    best = (0,0); best_n = 0
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            n = to_int(norm_func(a,b))
            if n <= 0: continue
            err = abs(n - m) / m
            if err < best_err:
                best_err = err; best = (a,b); best_n = n
    return best_n, best, float(best_err)

# ── Models ────────────────────────────────────────────────────────────────────
STRUCTURED = {
    "Z[omega]":   lambda a,b: norm_eisenstein(a,b),
    "Z[sqrt17]":  lambda a,b: norm_sqrt_d(a,b,17),
    "Z[sqrt5]":   lambda a,b: norm_sqrt_d(a,b,5),
    "Z[i]":       lambda a,b: norm_gaussian(a,b),
}

NULL_DS = {
    "Z[sqrt2]":   lambda a,b: norm_sqrt_d(a,b,2),
    "Z[sqrt3]":   lambda a,b: norm_sqrt_d(a,b,3),
    "Z[sqrt7]":   lambda a,b: norm_sqrt_d(a,b,7),
    "Z[sqrt11]":  lambda a,b: norm_sqrt_d(a,b,11),
    "Z[sqrt13]":  lambda a,b: norm_sqrt_d(a,b,13),
    "Z[sqrt19]":  lambda a,b: norm_sqrt_d(a,b,19),
    "Z[sqrt23]":  lambda a,b: norm_sqrt_d(a,b,23),
    "Z[sqrt29]":  lambda a,b: norm_sqrt_d(a,b,29),
    "Z[sqrt31]":  lambda a,b: norm_sqrt_d(a,b,31),
    "Z[sqrt37]":  lambda a,b: norm_sqrt_d(a,b,37),
}

def run(name, func, masses, R=150, verbose=False):
    results = {}
    if verbose: print(f"\n{name} (R={R}):")
    for f,m in masses.items():
        n,(a,b),err = find_best(m, func, R=R)
        results[f] = {'mass':float(m),'norm':n,'witness':(a,b),'error':float(err)}
        if verbose:
            print(f"  {f:>3}: m={m:>12.2f}  N={n:>10}  "
                  f"w=({a:4},{b:4})  err={err*100:.4f}%")
    return results

def stats(results):
    errs = [float(v['error']) for v in results.values()]
    return {
        'mean_log10': float(np.mean(np.log10(errs))),
        'median': float(np.median(errs)),
        'frac_2pct': float(sum(1 for e in errs if e<0.02)/len(errs)),
        'frac_5pct': float(sum(1 for e in errs if e<0.05)/len(errs)),
    }

# ── Run structured models ─────────────────────────────────────────────────────
print("="*70)
print("STRUCTURED MODELS (R=150)")
print("="*70)
struct_res = {}
struct_st  = {}
for name, func in STRUCTURED.items():
    res = run(name, func, MASSES, R=150, verbose=True)
    st  = stats(res)
    struct_res[name] = res
    struct_st[name]  = st

# ── Summary of structured models ─────────────────────────────────────────────
print()
print("="*70)
print("STRUCTURED MODEL SUMMARY")
print("="*70)
print(f"{'Model':<14} {'MLE':>8} {'median':>8} {'f<2%':>7} {'f<5%':>7}")
print("-"*45)
for name, st in struct_st.items():
    print(f"{name:<14} {st['mean_log10']:>8.3f} "
          f"{st['median']*100:>8.4f}% "
          f"{st['frac_2pct']:>7.2f} {st['frac_5pct']:>7.2f}")

# ── Run null models ───────────────────────────────────────────────────────────
print()
print("="*70)
print("NULL MODELS (R=100, no special arithmetic structure)")
print("="*70)
null_st = {}
for name, func in NULL_DS.items():
    res = run(name, func, MASSES, R=100, verbose=False)
    st  = stats(res)
    null_st[name] = st
    print(f"{name:<12}: MLE={st['mean_log10']:.3f}  "
          f"med={st['median']*100:.4f}%  "
          f"f<2%={st['frac_2pct']:.2f}  f<5%={st['frac_5pct']:.2f}")

# ── Null model averages ───────────────────────────────────────────────────────
null_mle_avg  = float(np.mean([s['mean_log10'] for s in null_st.values()]))
null_f2_avg   = float(np.mean([s['frac_2pct']  for s in null_st.values()]))
null_f5_avg   = float(np.mean([s['frac_5pct']  for s in null_st.values()]))
null_med_avg  = float(np.mean([s['median']      for s in null_st.values()]))

print(f"\nNull model average: MLE={null_mle_avg:.3f}  "
      f"med={null_med_avg*100:.4f}%  "
      f"f<2%={null_f2_avg:.2f}  f<5%={null_f5_avg:.2f}")

# ── Comparison ────────────────────────────────────────────────────────────────
print()
print("="*70)
print("COMPARISON: STRUCTURED vs NULL")
print("="*70)
print(f"{'Model':<14} {'MLE':>8} {'vs null':>8} {'f<2%':>7} {'vs null':>8} {'better?':>8}")
print("-"*60)
for name, st in struct_st.items():
    dm = st['mean_log10'] - null_mle_avg
    df = st['frac_2pct'] - null_f2_avg
    better = "YES" if dm < -0.3 or df > 0.15 else "marginal" if dm < 0 else "no"
    print(f"{name:<14} {st['mean_log10']:>8.3f} {dm:>+8.3f} "
          f"{st['frac_2pct']:>7.2f} {df:>+8.2f}  {better:>8}")

# ── Permutation test ──────────────────────────────────────────────────────────
print()
print("="*70)
print("PERMUTATION TEST (N=500 shuffles, R=80)")
print("="*70)
masses_list = list(MASSES.values())
keys_list   = list(MASSES.keys())

shuffle_f2  = {n: [] for n in STRUCTURED}
shuffle_mle = {n: [] for n in STRUCTURED}

for trial in range(500):
    random.seed(trial)
    sv = masses_list[:]
    random.shuffle(sv)
    sm = dict(zip(keys_list, sv))
    for name, func in STRUCTURED.items():
        res = run(name, func, sm, R=80, verbose=False)
        st  = stats(res)
        shuffle_f2[name].append(st['frac_2pct'])
        shuffle_mle[name].append(st['mean_log10'])

print(f"\n{'Model':<14} {'Obs f<2%':>9} {'Null mean':>10} "
      f"{'p(f)':>8} {'Obs MLE':>9} {'Null MLE':>9} {'p(MLE)':>8}")
print("-"*75)
for name in STRUCTURED:
    obs_f  = struct_st[name]['frac_2pct']
    obs_m  = struct_st[name]['mean_log10']
    nm_f   = float(np.mean(shuffle_f2[name]))
    nm_m   = float(np.mean(shuffle_mle[name]))
    pf     = float(sum(1 for x in shuffle_f2[name]  if x >= obs_f) / 500)
    pm     = float(sum(1 for x in shuffle_mle[name] if x <= obs_m) / 500)
    print(f"{name:<14} {obs_f:>9.2f} {nm_f:>10.2f} "
          f"{pf:>8.4f} {obs_m:>9.3f} {nm_m:>9.3f} {pm:>8.4f}")

# ── LaTeX table ───────────────────────────────────────────────────────────────
print()
print("="*70)
print("LATEX TABLE (best witnesses, two-field assignment)")
print("="*70)
print()
print(r"\begin{tabular}{lllrrr}")
print(r"\toprule")
print(r"Fermion & Sector & Ring & $m/m_e$ & Norm & Error \\")
print(r"\midrule")
for f in ["mu","tau"]:
    v = struct_res["Z[omega]"][f]
    a,b = v['witness']
    ring = r"$\mathbb{Z}[\omega]$"
    sector = "lepton"
    print(f"{f} & {sector} & {ring} & {v['mass']:.2f} & "
          f"{v['norm']} & ${v['error']*100:.4f}\\%$ \\\\")
for f in ["u","d","s","c","b","t"]:
    v = struct_res["Z[sqrt17]"][f]
    a,b = v['witness']
    ring = r"$\mathbb{Z}[\sqrt{17}]$"
    sector = "quark"
    print(f"{f} & {sector} & {ring} & {v['mass']:.2f} & "
          f"{v['norm']} & ${v['error']*100:.4f}\\%$ \\\\")
print(r"\bottomrule")
print(r"\end{tabular}")
