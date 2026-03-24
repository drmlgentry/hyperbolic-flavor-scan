"""
selberg_eigenvalue.py
=====================
Estimate the first Laplace eigenvalue lambda_1 for m003 and m006
using the Selberg zeta function evaluated on geodesic data.

Method
------
The Selberg zeta function for a compact hyperbolic 3-manifold is:

    Z(s) = prod_{gamma prim} prod_{m,n>=0}
           (1 - e^{i(m-n)phi(gamma)} * e^{-(s+m+n)*ell(gamma)})

Zeros of Z(s) on Re(s) = 1 are at s = 1 + i*r_k, where
lambda_k = s*(2-s) = 1 + r_k^2 are the Laplace eigenvalues.

We search for the smallest r > 0 where |Z(1+ir)| is minimised,
giving lambda_1 = 1 + r_1^2.

Practical truncation
--------------------
We use primitive geodesics up to MAX_LENGTH and truncate
the (m,n) product at m+n <= MN_MAX.  Larger MN_MAX gives
better accuracy but costs O(MN_MAX^2) per geodesic.

Reference
---------
Selberg (1956); Gangolli-Warner (1980); twist-angle paper DQ14199.

Usage
-----
    cd C:\dev\hyperbolic-flavor-scan
    python analysis\selberg_eigenvalue.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema
from scipy.optimize import minimize_scalar
import os, time

# ── Configuration ──────────────────────────────────────────────────
DATA_DIR   = r"C:\dev\hyperbolic-flavor-scan\data"
OUT_DIR    = r"C:\dev\hyperbolic-flavor-scan\analysis"

MANIFOLDS = {
    "m003": {"vol": 0.9814, "ref_lambda1": 2.48},
    "m006": {"vol": 2.0289, "ref_lambda1": 2.82},
}

MAX_LENGTH = 6.0   # use primitive geodesics up to this length
MN_MAX     = 4     # truncate m+n product at this value
R_MIN      = 0.5   # search for zero above this r
R_MAX      = 6.0   # search up to this r
R_POINTS   = 300   # scan resolution
MAX_PRIMS  = 2000  # use at most this many primitive geodesics


# ── Data loading ───────────────────────────────────────────────────
def load_geodesics(name: str) -> tuple[np.ndarray, np.ndarray]:
    """Load real lengths and twist angles from twist census CSV."""
    path = os.path.join(DATA_DIR, f"twist_census_len12_{name}.csv")
    df = pd.read_csv(path)
    cols = set(df.columns)

    # Real length
    if "real_length" in cols:
        ell = df["real_length"].values.astype(float)
    elif "abs_lambda" in cols:
        ell = 2.0 * np.log(np.abs(df["abs_lambda"].values.astype(float)))
    else:
        raise ValueError(f"No length column found. Columns: {list(cols)}")

    # Twist angle in radians
    if "phi_rad" in cols:
        phi = df["phi_rad"].values.astype(float)
    elif "phi_fold_deg" in cols:
        phi = np.radians(df["phi_fold_deg"].values.astype(float))
    elif "phi_deg" in cols:
        phi = np.radians(df["phi_deg"].values.astype(float))
    else:
        raise ValueError(f"No twist angle column found. Columns: {list(cols)}")

    # Filter: positive length, within MAX_LENGTH
    mask = np.isfinite(ell) & np.isfinite(phi) & (ell > 1e-8) & (ell <= MAX_LENGTH)
    return ell[mask], phi[mask]


def extract_primitives(ell: np.ndarray,
                       phi: np.ndarray,
                       tol: float = 1e-4) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract primitive geodesics.
    A geodesic is primitive if its length is not an integer multiple
    (within tol) of any shorter geodesic with consistent twist phase.

    For efficiency, work with unique (ell, phi) pairs sorted by ell.
    """
    # Deduplicate: keep unique (ell, phi_fold) pairs
    # phi_fold = min(|phi| mod pi, pi - |phi| mod pi) is the folded twist
    phi_fold = np.abs(phi) % np.pi
    phi_fold = np.minimum(phi_fold, np.pi - phi_fold)

    # Round to 4 decimal places and deduplicate
    ell_r   = np.round(ell, 4)
    phi_r   = np.round(phi_fold, 4)
    pairs   = np.unique(np.column_stack([ell_r, phi_r]), axis=0)

    # Sort by length
    order  = np.argsort(pairs[:, 0])
    pairs  = pairs[order]

    primitives = []
    for i, (l, p) in enumerate(pairs):
        is_prim = True
        for l0, p0 in pairs[:i]:
            if l0 < tol:
                continue
            k_f = l / l0
            k   = round(k_f)
            if k >= 2 and abs(k_f - k) < tol:
                # Check twist: phi ~ k * phi0 mod pi
                phase_diff = abs(p - k * p0) % np.pi
                phase_diff = min(phase_diff, np.pi - phase_diff)
                if phase_diff < tol:
                    is_prim = False
                    break
        if is_prim:
            primitives.append((l, p))

    prims = np.array(primitives)
    if len(prims) == 0:
        raise ValueError("No primitive geodesics found")
    return prims[:, 0], prims[:, 1]


# ── Selberg zeta ────────────────────────────────────────────────────
def log_selberg_zeta(r: float,
                     ell_prim: np.ndarray,
                     phi_prim: np.ndarray) -> complex:
    """
    Compute log|Z(1+ir)| using the truncated Selberg product.

    Z(s) = prod_{gamma prim} prod_{m,n>=0}
           (1 - e^{i(m-n)phi} * e^{-(s+m+n)*ell})

    with s = 1 + i*r, truncated at m+n <= MN_MAX.

    Returns Re(log Z) = sum of log|1 - ...| terms.
    """
    s  = 1.0 + 1j * r
    logZ = 0.0

    for ell, phi in zip(ell_prim, phi_prim):
        for m in range(MN_MAX + 1):
            for n in range(MN_MAX - m + 1):
                exp_arg = -(s + m + n) * ell + 1j * (m - n) * phi
                z = np.exp(exp_arg)
                denom = 1.0 - z
                if abs(denom) < 1e-14:
                    continue
                logZ += np.log(abs(denom))

    return logZ


def scan_zeta(ell_prim: np.ndarray,
              phi_prim: np.ndarray,
              r_vals: np.ndarray) -> np.ndarray:
    """Scan log|Z(1+ir)| over an array of r values."""
    log_amp = np.array([log_selberg_zeta(r, ell_prim, phi_prim) for r in r_vals])
    return log_amp


# ── Zero finding ────────────────────────────────────────────────────
def find_first_zero(ell_prim: np.ndarray,
                    phi_prim: np.ndarray,
                    r_vals: np.ndarray,
                    log_amp: np.ndarray) -> dict:
    """
    Find the first minimum of log|Z(1+ir)| above R_MIN.
    Refine with minimize_scalar.
    """
    # Smooth to reduce numerical noise
    smooth = gaussian_filter1d(log_amp, sigma=3)

    # Find local minima
    minima_idx = argrelextrema(smooth, np.less, order=5)[0]
    minima_idx = minima_idx[r_vals[minima_idx] > R_MIN]

    results = {}

    if len(minima_idx) == 0:
        results["r1_scan"]  = None
        results["lambda1"]  = None
        results["method"]   = "no minimum found"
        return results

    # First minimum
    i0    = minima_idx[0]
    r0    = r_vals[i0]
    r_lo  = max(R_MIN, r0 - 0.5)
    r_hi  = min(R_MAX, r0 + 0.5)

    results["r1_scan"]     = float(r0)
    results["logZ_at_r0"]  = float(log_amp[i0])

    # Refine
    try:
        res = minimize_scalar(
            lambda r: log_selberg_zeta(r, ell_prim, phi_prim),
            bounds=(r_lo, r_hi),
            method="bounded",
            options={"xatol": 1e-5, "maxiter": 100}
        )
        r1 = float(res.x)
        results["r1_refined"]  = r1
        results["lambda1"]     = 1.0 + r1**2
        results["logZ_min"]    = float(res.fun)
        results["method"]      = "minimize_scalar (bounded)"
    except Exception as e:
        results["r1_refined"] = r0
        results["lambda1"]    = 1.0 + r0**2
        results["method"]     = f"scan only (refine failed: {e})"

    # List all minima found
    results["all_minima_r"] = [float(r_vals[i]) for i in minima_idx[:5]]
    results["all_minima_lambda"] = [1.0 + r_vals[i]**2 for i in minima_idx[:5]]

    return results


# ── Plotting ────────────────────────────────────────────────────────
def plot_zeta(name: str,
              r_vals: np.ndarray,
              log_amp: np.ndarray,
              fit: dict):
    """Plot log|Z(1+ir)| with zero locations marked."""
    fig, ax = plt.subplots(figsize=(10, 5))

    smooth = gaussian_filter1d(log_amp, sigma=3)
    ax.plot(r_vals, log_amp, "b-",  lw=0.8, alpha=0.5, label="log|Z(1+ir)|")
    ax.plot(r_vals, smooth,  "b-",  lw=2,   label="smoothed")

    # Mark zeros
    if fit.get("r1_refined"):
        r1  = fit["r1_refined"]
        lam = fit["lambda1"]
        ax.axvline(r1, color="red", lw=2, ls="--",
                   label=fr"$r_1 = {r1:.3f}$, $\lambda_1 = {lam:.3f}$")

    # Reference
    ref = MANIFOLDS[name]["ref_lambda1"]
    r_ref = np.sqrt(max(ref - 1, 0))
    ax.axvline(r_ref, color="green", lw=1.5, ls=":",
               label=fr"Reference $r_1 = {r_ref:.3f}$ ($\lambda_1 = {ref}$)")

    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xlabel(r"$r$", fontsize=12)
    ax.set_ylabel(r"$\log|Z(1+ir)|$", fontsize=12)
    ax.set_title(f"{name} — Selberg zeta scan  "
                 f"(L_max={MAX_LENGTH}, m+n≤{MN_MAX}, "
                 f"{len(fit.get('all_minima_r', []))} minima)",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    outpath = os.path.join(OUT_DIR, f"selberg_zeta_{name}.png")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {outpath}")


# ── Main ────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 62)
    print("Selberg zeta eigenvalue estimation for m003 and m006")
    print(f"MAX_LENGTH={MAX_LENGTH}  MN_MAX={MN_MAX}  "
          f"R=[{R_MIN},{R_MAX}]  points={R_POINTS}")
    print("=" * 62)

    r_vals = np.linspace(R_MIN, R_MAX, R_POINTS)

    for name in ["m003", "m006"]:
        print(f"\n{'─'*40}")
        print(f"Manifold: {name}  (vol={MANIFOLDS[name]['vol']:.4f})")

        # Load
        t0 = time.time()
        ell_all, phi_all = load_geodesics(name)
        print(f"  Loaded {len(ell_all)} geodesics (ell <= {MAX_LENGTH})")

        # Extract primitives
        ell_p, phi_p = extract_primitives(ell_all, phi_all)
        # Limit for speed
        if len(ell_p) > MAX_PRIMS:
            # Keep shortest (most important for convergence)
            order = np.argsort(ell_p)
            ell_p = ell_p[order[:MAX_PRIMS]]
            phi_p = phi_p[order[:MAX_PRIMS]]
        print(f"  Primitive geodesics used: {len(ell_p)}")
        print(f"  Systole: {ell_p.min():.6f}")

        # Scan
        print(f"  Scanning r in [{R_MIN}, {R_MAX}] ({R_POINTS} points)...")
        log_amp = scan_zeta(ell_p, phi_p, r_vals)
        t1 = time.time()
        print(f"  Scan time: {t1-t0:.1f}s")

        # Find zero
        fit = find_first_zero(ell_p, phi_p, r_vals, log_amp)

        # Report
        print(f"\n  Results:")
        print(f"    Scan minimum at r = {fit.get('r1_scan', 'N/A')}")
        if fit.get("r1_refined"):
            print(f"    Refined minimum: r_1 = {fit['r1_refined']:.5f}")
            print(f"    lambda_1 = {fit['lambda1']:.5f}")
            print(f"    log|Z| at minimum = {fit.get('logZ_min', 'N/A'):.4f}")
        print(f"    Method: {fit.get('method', 'N/A')}")
        print(f"    All minima r: {[f'{x:.3f}' for x in fit.get('all_minima_r', [])]}")
        print(f"    All minima lambda: {[f'{x:.3f}' for x in fit.get('all_minima_lambda', [])]}")
        ref = MANIFOLDS[name]["ref_lambda1"]
        print(f"    Reference (twist paper): lambda_1 ~ {ref}")

        # Plot
        plot_zeta(name, r_vals, log_amp, fit)

    print(f"\n{'='*62}")
    print("IMPORTANT CAVEATS")
    print("  1. Truncation at L_max and MN_MAX introduces systematic error.")
    print("  2. Results should be compared to the twist paper values (DQ14199).")
    print("  3. The minimum of log|Z| is not necessarily a zero; a true zero")
    print("     requires log|Z| -> -inf.  The truncated product gives a finite")
    print("     minimum whose location approximates the true zero.")
    print("  4. For rigorous values, use the finite element method or")
    print("     extend the geodesic sum to L_max > 20.")


if __name__ == "__main__":
    main()
