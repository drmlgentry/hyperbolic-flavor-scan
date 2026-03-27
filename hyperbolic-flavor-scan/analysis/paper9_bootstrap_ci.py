import pandas as pd, numpy as np
from scipy.optimize import curve_fit

df = pd.read_csv(r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m006.csv")

print("Bootstrap CIs on exponential decay rates c_k (word length basis)")
print("="*60)

N_BOOTSTRAP = 2000

for k in range(5):
    sub = df[df.h1_class == k]
    lengths = sorted(sub["length"].unique())
    floors  = []
    for L in lengths:
        window = sub[sub.length <= L]
        f = window.phi_fold_deg.min()
        if f > 1e-6:
            floors.append((L, f))
    
    if len(floors) < 4:
        print(f"Class {k}: insufficient data"); continue
    
    Ls = np.array([x[0] for x in floors])
    Fs = np.array([x[1] for x in floors])
    
    def exp_f(x, A, c): return A * np.exp(-c * x)
    
    try:
        popt, _ = curve_fit(exp_f, Ls, Fs, p0=[50, 0.7], maxfev=5000)
        c_hat = popt[1]
    except:
        print(f"Class {k}: fit failed"); continue
    
    c_boots = []
    rng = np.random.default_rng(42 + k)
    for _ in range(N_BOOTSTRAP):
        idx = rng.integers(0, len(Ls), len(Ls))
        try:
            pb, _ = curve_fit(exp_f, Ls[idx], Fs[idx], p0=[50, 0.7], maxfev=2000)
            if 0 < pb[1] < 5: c_boots.append(pb[1])
        except: pass
    
    if len(c_boots) < 100:
        print(f"Class {k}: bootstrap unstable ({len(c_boots)} converged)"); continue
    
    c_boots = np.array(c_boots)
    lo, hi  = np.percentile(c_boots, [2.5, 97.5])
    print(f"Class {k}: c={c_hat:.3f}  95% CI [{lo:.3f}, {hi:.3f}]  "
          f"({len(c_boots)}/{N_BOOTSTRAP} converged)")

print("\nPairing check: c_1 CI and c_4 CI should overlap; same for c_2/c_3.")
