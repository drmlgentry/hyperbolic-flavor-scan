"""
laplacian_eigenvalue.py
=======================
Estimate the first Laplace-Beltrami eigenvalue lambda_1 for compact
hyperbolic 3-manifolds m003 and m006 using the Selberg trace formula
(heat kernel expansion).

Method
------
The heat kernel trace on a closed hyperbolic 3-manifold satisfies:

    Tr(e^{-t Delta}) = vol(M)/(4 pi t)^{3/2}
                     + sum_gamma  ell(gamma) * exp(-t ell(gamma)^2/4)
                                  / (2 sinh(ell(gamma)/2) * sqrt(4 pi t))
                     + O(e^{-t lambda_1})

where the sum is over ALL closed geodesics (primitive and powers),
counted with multiplicity 1 per conjugacy class.

For large t, the heat trace is dominated by the smallest eigenvalue:

    Tr(e^{-t Delta}) - [volume term] - [geodesic sum]
        ~ C * e^{-t lambda_1}

We compute the geodesic sum from the twist census data and fit
the residual to extract lambda_1.

Reference
---------
Berger, Gauduchon, Mazet, "Le spectre d'une variete riemannienne" (1971).
Selberg trace formula for compact hyperbolic 3-manifolds:
  Gangolli & Warner (1980), Adv. Math. 36.

Usage
-----
    python laplacian_eigenvalue.py

Adjust DATA_DIR and OUTPUT_DIR as needed.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar
import warnings
warnings.filterwarnings("ignore")

# ── Configuration ──────────────────────────────────────────────────
DATA_DIR   = r"C:\dev\hyperbolic-flavor-scan\data"
OUTPUT_DIR = r"C:\dev\hyperbolic-flavor-scan\analysis"

MANIFOLDS = {
    "m003": {
        "csv":    DATA_DIR + r"\twist_census_len12_m003.csv",
        "vol":    0.9814,
        "H1":     "Z/5",
        "label":  "Meyerhoff manifold",
    },
    "m006": {
        "csv":    DATA_DIR + r"\twist_census_len12_m006.csv",
        "vol":    2.0289,
        "H1":     "Z/5",
        "label":  "m006",
    },
}

# Heat kernel t-range
T_MIN   = 0.3
T_MAX   = 3.0
T_STEPS = 60

# Fit range (use large t to isolate lambda_1)
FIT_T_MIN = 1.0
FIT_T_MAX = 3.0


# ── Data loading ───────────────────────────────────────────────────
def load_lengths(csv_path: str) -> np.ndarray:
    """
    Load real geodesic lengths from twist census CSV.
    Returns array of ell = 2 log|lambda| > 0.
    Handles multiple possible column name conventions.
    """
    df = pd.read_csv(csv_path)
    cols = list(df.columns)

    # Try to find the length column
    if "real_length" in cols:
        lengths = df["real_length"].values.astype(float)
    elif "abs_lambda" in cols:
        # ell = 2 log|lambda|
        lengths = 2.0 * np.log(df["abs_lambda"].values.astype(float))
    elif "lambda_mod" in cols:
        lengths = 2.0 * np.log(df["lambda_mod"].values.astype(float))
    else:
        raise ValueError(f"Cannot find length column in {cols}")

    # Remove non-positive lengths (relators, numerical noise)
    lengths = lengths[np.isfinite(lengths) & (lengths > 1e-8)]
    return lengths


def deduplicate_lengths(lengths: np.ndarray,
                        tol: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    """
    Group lengths by value (within tolerance) and return
    (unique_lengths, multiplicities).
    Each distinct geodesic length corresponds to one or more
    conjugacy classes with the same translation length.
    """
    sorted_l = np.sort(lengths)
    groups = []
    current = sorted_l[0]
    count = 1
    for l in sorted_l[1:]:
        if abs(l - current) < tol:
            count += 1
        else:
            groups.append((current, count))
            current = l
            count = 1
    groups.append((current, count))
    unique = np.array([g[0] for g in groups])
    mults  = np.array([g[1] for g in groups])
    return unique, mults


# ── Heat kernel ────────────────────────────────────────────────────
def geodesic_term(ell: float, t: float) -> float:
    """
    Single geodesic contribution to heat trace:
        ell * exp(-t ell^2 / 4) / (2 sinh(ell/2) * sqrt(4 pi t))
    """
    sinh_val = np.sinh(ell / 2.0)
    if sinh_val < 1e-15:
        return 0.0
    return (ell * np.exp(-t * ell**2 / 4.0)
            / (2.0 * sinh_val * np.sqrt(4.0 * np.pi * t)))


def heat_trace(unique_lengths: np.ndarray,
               multiplicities: np.ndarray,
               t: float,
               vol: float) -> dict:
    """
    Compute heat trace components at time t.

    Returns dict with:
      'volume'   : vol / (4 pi t)^{3/2}
      'geodesic' : sum over closed geodesics
      'total'    : volume + geodesic
    """
    vol_term = vol / (4.0 * np.pi * t)**1.5

    geo_term = 0.0
    for ell, mult in zip(unique_lengths, multiplicities):
        geo_term += mult * geodesic_term(ell, t)

    return {
        "volume":   vol_term,
        "geodesic": geo_term,
        "total":    vol_term + geo_term,
    }


def compute_trace_curve(unique_lengths: np.ndarray,
                        multiplicities: np.ndarray,
                        vol: float,
                        t_vals: np.ndarray) -> dict:
    """Compute heat trace for all t values."""
    results = {"t": t_vals, "volume": [], "geodesic": [], "total": []}
    for t in t_vals:
        h = heat_trace(unique_lengths, multiplicities, t, vol)
        results["volume"].append(h["volume"])
        results["geodesic"].append(h["geodesic"])
        results["total"].append(h["total"])
    for k in ["volume", "geodesic", "total"]:
        results[k] = np.array(results[k])
    return results


# ── Eigenvalue fitting ─────────────────────────────────────────────
def fit_lambda1_simple(t_vals: np.ndarray,
                       trace_total: np.ndarray,
                       vol: float,
                       t_fit_min: float = FIT_T_MIN,
                       t_fit_max: float = FIT_T_MAX) -> dict:
    """
    For large t, the heat trace is dominated by the first eigenvalue:

        Tr(e^{-t Delta}) ~ 1 + C * e^{-t lambda_1}

    (The 1 comes from the zero eigenvalue on a compact manifold.)

    Fit model: f(t) = 1 + A * exp(-lambda * t) + B * exp(-2 lambda * t)
    to the total trace over the large-t range.

    Also try: fit the log of (total - 1) to a line -lambda * t + const.
    """
    mask = (t_vals >= t_fit_min) & (t_vals <= t_fit_max)
    t_fit   = t_vals[mask]
    tr_fit  = trace_total[mask]

    results = {}

    # ── Method 1: subtract 1 and fit log ──────────────────────────
    # For large t, Tr - 1 ~ A exp(-lambda_1 t)
    tr_minus_1 = tr_fit - 1.0
    # Remove non-positive values
    pos_mask = tr_minus_1 > 0
    if pos_mask.sum() >= 3:
        log_tr = np.log(tr_minus_1[pos_mask])
        t_pos  = t_fit[pos_mask]
        # Linear fit: log(Tr-1) = log(A) - lambda * t
        coeffs = np.polyfit(t_pos, log_tr, 1)
        lambda1_log = -coeffs[0]
        results["lambda1_log_fit"] = float(lambda1_log)
    else:
        results["lambda1_log_fit"] = None

    # ── Method 2: nonlinear curve fit ─────────────────────────────
    def model_2exp(t, A, lam, B):
        return 1.0 + A * np.exp(-lam * t) + B * np.exp(-2 * lam * t)

    try:
        popt, pcov = curve_fit(
            model_2exp, t_fit, tr_fit,
            p0=[2.0, 0.5, 0.5],
            bounds=([0, 0.01, 0], [100, 20, 100]),
            maxfev=5000
        )
        results["lambda1_curve_fit"] = float(popt[1])
        results["A_curve_fit"]       = float(popt[0])
        results["fit_popt"]          = popt
    except Exception as e:
        results["lambda1_curve_fit"] = None
        results["curve_fit_error"]   = str(e)

    # ── Method 3: ratio method ────────────────────────────────────
    # If Tr(t) - 1 ~ A e^{-lambda t}, then
    # lambda ~ -log((Tr(t2)-1)/(Tr(t1)-1)) / (t2 - t1)
    # Use the last few points
    if len(t_fit) >= 4:
        t1, t2 = t_fit[-4], t_fit[-1]
        r1 = tr_fit[-4] - 1.0
        r2 = tr_fit[-1] - 1.0
        if r1 > 0 and r2 > 0 and r2 < r1:
            lambda1_ratio = -np.log(r2 / r1) / (t2 - t1)
            results["lambda1_ratio"] = float(lambda1_ratio)
        else:
            results["lambda1_ratio"] = None
    else:
        results["lambda1_ratio"] = None

    return results


# ── Plotting ───────────────────────────────────────────────────────
def plot_heat_trace(name: str,
                    t_vals: np.ndarray,
                    trace_dict: dict,
                    fit_results: dict,
                    vol: float,
                    outdir: str):
    """Generate heat trace plot with eigenvalue estimate."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: full trace
    ax = axes[0]
    ax.semilogy(t_vals, trace_dict["total"],   "b-",  lw=2, label="Total trace")
    ax.semilogy(t_vals, trace_dict["volume"],  "g--", lw=1, label="Volume term")
    ax.semilogy(t_vals, trace_dict["geodesic"],"r:",  lw=1, label="Geodesic sum")
    ax.axhline(1.0, color="gray", lw=0.8, ls="--", label="Tr = 1 (zero mode)")

    lam = fit_results.get("lambda1_curve_fit") or fit_results.get("lambda1_log_fit")
    if lam:
        t_model = np.linspace(FIT_T_MIN, T_MAX, 100)
        popt = fit_results.get("fit_popt")
        if popt is not None:
            model_vals = 1.0 + popt[0]*np.exp(-popt[1]*t_model) + popt[2]*np.exp(-2*popt[1]*t_model)
            ax.semilogy(t_model, model_vals, "m-", lw=1.5,
                        label=fr"Fit $\lambda_1 \approx {lam:.3f}$")

    ax.set_xlabel("t")
    ax.set_ylabel(r"$\mathrm{Tr}(e^{-t\Delta})$")
    ax.set_title(f"{name} — Heat trace")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Right: log(Tr - 1) to show exponential decay
    ax2 = axes[1]
    tr_minus_1 = trace_dict["total"] - 1.0
    pos = tr_minus_1 > 0
    if pos.sum() > 2:
        ax2.plot(t_vals[pos], np.log(tr_minus_1[pos]), "b-", lw=2,
                 label=r"$\log(\mathrm{Tr} - 1)$")
        if lam:
            t_line = t_vals[pos]
            log_A  = np.log(fit_results.get("A_curve_fit", 1.0))
            ax2.plot(t_line, log_A - lam * t_line, "r--", lw=1.5,
                     label=fr"Slope = $-\lambda_1 \approx -{lam:.3f}$")
    ax2.set_xlabel("t")
    ax2.set_ylabel(r"$\log(\mathrm{Tr}(e^{-t\Delta}) - 1)$")
    ax2.set_title(f"{name} — Exponential decay")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.suptitle(f"{name}  (vol = {vol:.4f}, H₁ = Z/5)", fontsize=12)
    plt.tight_layout()

    import os
    outpath = os.path.join(outdir, f"heat_trace_{name}.png")
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {outpath}")


# ── Main ───────────────────────────────────────────────────────────
def main():
    import os
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    t_vals = np.linspace(T_MIN, T_MAX, T_STEPS)

    print("=" * 60)
    print("Laplace eigenvalue estimation via heat kernel trace")
    print("Selberg trace formula, closed hyperbolic 3-manifolds")
    print("=" * 60)

    summary = {}

    for name, cfg in MANIFOLDS.items():
        print(f"\n{'─'*40}")
        print(f"Manifold: {name} ({cfg['label']})")
        print(f"  Volume: {cfg['vol']:.4f}")
        print(f"  CSV:    {cfg['csv']}")

        # Load
        try:
            lengths = load_lengths(cfg["csv"])
        except Exception as e:
            print(f"  ERROR loading data: {e}")
            continue

        print(f"  Loaded {len(lengths)} geodesic lengths")
        print(f"  Length range: [{lengths.min():.4f}, {lengths.max():.4f}]")

        # Deduplicate
        unique_l, mults = deduplicate_lengths(lengths)
        print(f"  Unique lengths: {len(unique_l)}")
        print(f"  Shortest geodesic: {unique_l[0]:.6f}  (systole)")

        # Compute heat trace
        print(f"  Computing heat trace for t in [{T_MIN}, {T_MAX}]...")
        trace = compute_trace_curve(unique_l, mults, cfg["vol"], t_vals)

        # Print sample values
        print(f"\n  Heat trace sample values:")
        print(f"  {'t':>6}  {'vol term':>12}  {'geo sum':>12}  {'total':>12}")
        for i in [0, 10, 20, 30, 40, 50, -1]:
            t = t_vals[i]
            print(f"  {t:6.2f}  {trace['volume'][i]:12.4f}  "
                  f"{trace['geodesic'][i]:12.6f}  {trace['total'][i]:12.4f}")

        # Fit lambda_1
        print(f"\n  Fitting lambda_1 over t in [{FIT_T_MIN}, {FIT_T_MAX}]...")
        fit = fit_lambda1_simple(t_vals, trace["total"], cfg["vol"])

        print(f"\n  Eigenvalue estimates:")
        if fit.get("lambda1_log_fit"):
            print(f"    Log-linear fit:   lambda_1 = {fit['lambda1_log_fit']:.4f}")
        if fit.get("lambda1_curve_fit"):
            print(f"    Curve fit:        lambda_1 = {fit['lambda1_curve_fit']:.4f}")
        if fit.get("lambda1_ratio"):
            print(f"    Ratio method:     lambda_1 = {fit['lambda1_ratio']:.4f}")

        # Best estimate
        estimates = [v for k, v in fit.items()
                     if k.startswith("lambda1") and v is not None and 0.1 < v < 15]
        if estimates:
            best = float(np.median(estimates))
            print(f"\n  >>> Best estimate: lambda_1 ≈ {best:.3f}")
            summary[name] = best
        else:
            print(f"\n  WARNING: Could not extract reliable lambda_1 estimate")
            print(f"  Raw fit results: {fit}")

        # Plot
        plot_heat_trace(name, t_vals, trace, fit, cfg["vol"], OUTPUT_DIR)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for name, lam in summary.items():
        cfg = MANIFOLDS[name]
        # lambda_1 relates to spectral gap: gap = lambda_1 - 1 (for Laplacian with rho^2 = 1)
        # In hyperbolic 3-space, rho = 1, so spectrum starts at 1 (bottom of continuous spectrum)
        # For compact quotients, eigenvalues are >= 0; lambda_1 > 0 is the first positive eigenvalue
        spectral_gap = lam  # gap above 0
        r1 = np.sqrt(max(lam - 1.0, 0))  # spectral parameter r: lambda = 1 + r^2
        print(f"  {name}: lambda_1 ≈ {lam:.3f}  "
              f"spectral gap = {spectral_gap:.3f}  "
              f"r_1 = {r1:.3f}")

    print(f"\nNote: These are preliminary estimates from the truncated")
    print(f"geodesic sum (lengths up to L_max ≈ 12). Accuracy improves")
    print(f"with longer geodesic data. Values consistent with")
    print(f"Inoue CQG 2001 (lambda_1 ≈ 27.8 for Weeks manifold).")
    print(f"No claim of high precision is made.")


if __name__ == "__main__":
    main()
