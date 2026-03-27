"""
selberg_app.py
==============
Streamlit UI for Selberg zeta eigenvalue computation
with real-time progress, checkpointing, and results visualization.

Run:
    streamlit run selberg_app.py
"""

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.ndimage import gaussian_filter1d
from scipy.signal import argrelextrema
import time, os, json, hashlib, io
from datetime import datetime
from pathlib import Path

# ── Page config ────────────────────────────────────────────────────
st.set_page_config(
    page_title="Selberg Zeta · Eigenvalue Engine",
    page_icon="λ",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Stylesheet ─────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;700&family=Libre+Baskerville:ital,wght@0,400;0,700;1,400&display=swap');

:root {
  --bg:      #0a0e17;
  --surface: #111827;
  --border:  #1f2d3d;
  --accent:  #4fc3f7;
  --accent2: #81c784;
  --warn:    #ffb74d;
  --text:    #cdd5e0;
  --dim:     #6b7280;
  --red:     #ef5350;
}

html, body, [class*="css"] {
  background-color: var(--bg) !important;
  color: var(--text) !important;
  font-family: 'Libre Baskerville', Georgia, serif;
}

.stApp { background-color: var(--bg) !important; }

h1, h2, h3 {
  font-family: 'Libre Baskerville', serif !important;
  color: var(--accent) !important;
  letter-spacing: 0.02em;
}

.mono { font-family: 'JetBrains Mono', monospace; }

.metric-card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-left: 3px solid var(--accent);
  border-radius: 6px;
  padding: 14px 18px;
  margin: 6px 0;
}

.eigenvalue-card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 16px;
  text-align: center;
  font-family: 'JetBrains Mono', monospace;
}

.eigenvalue-card .lam {
  font-size: 2em;
  color: var(--accent);
  font-weight: 700;
}

.eigenvalue-card .sub {
  color: var(--dim);
  font-size: 0.85em;
  margin-top: 4px;
}

.log-box {
  background: #060a10;
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 12px 14px;
  font-family: 'JetBrains Mono', monospace;
  font-size: 0.78em;
  color: #8bc34a;
  max-height: 300px;
  overflow-y: auto;
  white-space: pre-wrap;
}

.progress-step {
  display: flex;
  align-items: center;
  gap: 10px;
  padding: 6px 0;
  font-family: 'JetBrains Mono', monospace;
  font-size: 0.85em;
}

.step-done  { color: var(--accent2); }
.step-run   { color: var(--warn);    }
.step-wait  { color: var(--dim);     }

.checkpoint-badge {
  display: inline-block;
  background: #1a2a1a;
  color: var(--accent2);
  border: 1px solid var(--accent2);
  border-radius: 4px;
  padding: 2px 8px;
  font-family: 'JetBrains Mono', monospace;
  font-size: 0.75em;
}

.warning-box {
  background: #1a1200;
  border: 1px solid var(--warn);
  border-radius: 6px;
  padding: 10px 14px;
  font-size: 0.85em;
  color: var(--warn);
}
</style>
""", unsafe_allow_html=True)


# ── Defaults ───────────────────────────────────────────────────────
DEFAULT_CSV = {
    "m003": r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m003.csv",
    "m006": r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m006.csv",
}
DEFAULT_CHECKPOINT_DIR = r"C:\dev\hyperbolic-flavor-scan\analysis\checkpoints"
VOLUMES = {"m003": 0.9814, "m006": 2.0289}
REF_LAMBDA = {"m003": 2.48, "m006": 2.82}
PHI = (1 + np.sqrt(5)) / 2


# ── Checkpoint system ──────────────────────────────────────────────
def checkpoint_path(ckpt_dir: str, name: str, params: dict) -> Path:
    key = json.dumps(params, sort_keys=True)
    h = hashlib.md5(key.encode()).hexdigest()[:8]
    return Path(ckpt_dir) / f"{name}_{h}.json"

def save_checkpoint(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)

def load_checkpoint(path: Path) -> dict | None:
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


# ── Core computation ───────────────────────────────────────────────
@st.cache_data(show_spinner=False)
def load_geodesics(csv_path: str, max_length: float):
    df = pd.read_csv(csv_path)
    cols = set(df.columns)
    if "real_length" in cols:
        ell = df["real_length"].values.astype(float)
    elif "abs_lambda" in cols:
        ell = 2.0 * np.log(np.abs(df["abs_lambda"].values.astype(float)))
    else:
        raise ValueError(f"No length column. Columns: {list(cols)}")
    if "phi_rad" in cols:
        phi = df["phi_rad"].values.astype(float)
    elif "phi_fold_deg" in cols:
        phi = np.radians(df["phi_fold_deg"].values.astype(float))
    elif "phi_deg" in cols:
        phi = np.radians(df["phi_deg"].values.astype(float))
    else:
        raise ValueError(f"No twist column. Columns: {list(cols)}")
    mask = np.isfinite(ell) & np.isfinite(phi) & (ell > 1e-8) & (ell <= max_length)
    return ell[mask], phi[mask]


@st.cache_data(show_spinner=False)
def extract_primitives(ell: np.ndarray, phi: np.ndarray,
                       max_prims: int, tol: float = 1e-4):
    phi_fold = np.abs(phi) % np.pi
    phi_fold = np.minimum(phi_fold, np.pi - phi_fold)
    ell_r = np.round(ell, 4)
    phi_r = np.round(phi_fold, 4)
    pairs = np.unique(np.column_stack([ell_r, phi_r]), axis=0)
    order = np.argsort(pairs[:, 0])
    pairs = pairs[order]

    primitives = []
    for i, (l, p) in enumerate(pairs):
        is_prim = True
        for l0, p0 in pairs[:i]:
            if l0 < tol: continue
            k_f = l / l0
            k = round(k_f)
            if k >= 2 and abs(k_f - k) < tol:
                pd_ = abs(p - k * p0) % np.pi
                pd_ = min(pd_, np.pi - pd_)
                if pd_ < tol:
                    is_prim = False
                    break
        if is_prim:
            primitives.append((l, p))

    prims = np.array(primitives) if primitives else np.zeros((0, 2))
    if len(prims) > max_prims:
        prims = prims[:max_prims]
    return prims[:, 0], prims[:, 1]


def log_zeta(r: float, ell_p: np.ndarray, phi_p: np.ndarray,
             mn_max: int) -> float:
    s = 1.0 + 1j * r
    logZ = 0.0
    for ell, phi in zip(ell_p, phi_p):
        for m in range(mn_max + 1):
            for n in range(mn_max - m + 1):
                exp_arg = -(s + m + n) * ell + 1j * (m - n) * phi
                z = np.exp(exp_arg)
                d = 1.0 - z
                if abs(d) > 1e-14:
                    logZ += np.log(abs(d))
    return logZ


def run_scan(ell_p, phi_p, mn_max, r_vals,
             progress_cb=None, log_cb=None) -> np.ndarray:
    """Scan log|Z(1+ir)| over r_vals with progress callbacks."""
    log_amp = np.full(len(r_vals), np.nan)
    n = len(r_vals)
    t0 = time.time()
    for i, r in enumerate(r_vals):
        log_amp[i] = log_zeta(r, ell_p, phi_p, mn_max)
        if progress_cb and i % max(1, n // 200) == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (n - i - 1)
            progress_cb(i / n, elapsed, eta)
        if log_cb and i % max(1, n // 50) == 0:
            log_cb(f"  r={r:.3f}  log|Z|={log_amp[i]:.4f}  "
                   f"({i}/{n}  elapsed {elapsed:.0f}s  ETA {eta:.0f}s)")
    return log_amp


def find_minima(r_vals, log_amp, ell_p, phi_p, mn_max,
                n_eigenvalues=5) -> list[dict]:
    smooth = gaussian_filter1d(log_amp, sigma=4)
    minima_idx = argrelextrema(smooth, np.less, order=8)[0]
    results = []
    seen_r = []
    for idx in minima_idx:
        r0 = r_vals[idx]
        # Skip if too close to a previous minimum
        if any(abs(r0 - rs) < 0.3 for rs in seen_r):
            continue
        r_lo = max(r_vals[0], r0 - 0.4)
        r_hi = min(r_vals[-1], r0 + 0.4)
        try:
            res = minimize_scalar(
                lambda r: log_zeta(r, ell_p, phi_p, mn_max),
                bounds=(r_lo, r_hi), method="bounded",
                options={"xatol": 1e-5, "maxiter": 80}
            )
            r1 = float(res.x)
            lam = 1.0 + r1 ** 2
            logz = float(res.fun)
            seen_r.append(r1)
            results.append({"r": r1, "lambda": lam, "logZ_min": logz})
        except Exception:
            pass
        if len(results) >= n_eigenvalues:
            break
    return sorted(results, key=lambda x: x["lambda"])


def phi_lattice_check(lam: float) -> dict:
    """Check if lambda is near a phi-lattice point."""
    me = 5.11e-4  # GeV
    # Convert lambda to an effective mass scale?
    # Actually check lambda itself against the lattice
    log_lam = np.log(lam) / np.log(PHI)
    q_star = round(4 * log_lam)
    residual = abs(log_lam - q_star / 4)
    predicted = PHI ** (q_star / 4)
    return {"q": q_star, "residual": residual,
            "predicted": predicted, "pct_err": abs(lam - predicted) / lam * 100}


# ── Plot ───────────────────────────────────────────────────────────
def make_plot(name, r_vals, log_amp, minima, ref_lam):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                              facecolor="#0a0e17")
    for ax in axes:
        ax.set_facecolor("#0a0e17")
        ax.tick_params(colors="#6b7280")
        for spine in ax.spines.values():
            spine.set_color("#1f2d3d")

    smooth = gaussian_filter1d(log_amp, sigma=4)

    # Left: log|Z| scan
    ax = axes[0]
    ax.plot(r_vals, log_amp, color="#1f3a5c", lw=0.7, alpha=0.6)
    ax.plot(r_vals, smooth,  color="#4fc3f7", lw=1.8, label="log|Z(1+ir)|")
    ax.axhline(0, color="#1f2d3d", lw=0.8)

    ref_r = np.sqrt(max(ref_lam - 1, 0))
    ax.axvline(ref_r, color="#ffb74d", lw=1.2, ls=":",
               label=f"Reference r={ref_r:.3f} (λ={ref_lam})")

    colors = ["#ef5350", "#81c784", "#ce93d8", "#ffcc02", "#4dd0e1"]
    for i, m in enumerate(minima[:5]):
        c = colors[i % len(colors)]
        ax.axvline(m["r"], color=c, lw=1.5, ls="--", alpha=0.9)
        ax.annotate(f"λ={m['lambda']:.3f}",
                    xy=(m["r"], smooth[np.argmin(np.abs(r_vals - m["r"]))]),
                    xytext=(m["r"] + 0.05, smooth.max() * 0.7 - i * smooth.std() * 0.3),
                    color=c, fontsize=7, fontfamily="monospace",
                    arrowprops=dict(arrowstyle="-", color=c, lw=0.8))

    ax.set_xlabel("r", color="#cdd5e0")
    ax.set_ylabel("log|Z(1+ir)|", color="#cdd5e0")
    ax.set_title(f"{name}  —  Selberg zeta scan",
                 color="#4fc3f7", fontsize=11)
    ax.legend(fontsize=8, facecolor="#111827",
              labelcolor="#cdd5e0", edgecolor="#1f2d3d")
    ax.grid(True, color="#0f1a28", lw=0.5)

    # Right: eigenvalue spectrum
    ax2 = axes[1]
    if minima:
        lams = [m["lambda"] for m in minima]
        bars = ax2.barh(range(len(lams)), lams,
                        color=[colors[i % len(colors)] for i in range(len(lams))],
                        height=0.5, alpha=0.85)
        ax2.set_yticks(range(len(lams)))
        ax2.set_yticklabels([f"λ_{i+1}" for i in range(len(lams))],
                             color="#cdd5e0", fontfamily="monospace")
        ax2.set_xlabel("λ = 1 + r²", color="#cdd5e0")
        ax2.set_title("Low-lying eigenvalues", color="#4fc3f7", fontsize=11)
        ax2.axvline(ref_lam, color="#ffb74d", lw=1.2, ls=":",
                    label=f"Reference λ₁={ref_lam}")
        ax2.axvline(0.75, color="#6b7280", lw=0.8, ls="--",
                    label="Selberg threshold 3/4")
        for bar, lam in zip(bars, lams):
            ax2.text(lam + 0.05, bar.get_y() + bar.get_height() / 2,
                     f"{lam:.4f}", va="center", color="#cdd5e0",
                     fontsize=8, fontfamily="monospace")
        ax2.legend(fontsize=8, facecolor="#111827",
                   labelcolor="#cdd5e0", edgecolor="#1f2d3d")
        ax2.grid(True, color="#0f1a28", lw=0.5, axis="x")

    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=160,
                bbox_inches="tight", facecolor="#0a0e17")
    buf.seek(0)
    plt.close()
    return buf


# ── Sidebar ────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## λ Selberg Eigenvalue Engine")
    st.markdown("---")

    manifold = st.selectbox("Manifold", ["m003 (Meyerhoff)", "m006"])
    name = manifold.split()[0]

    csv_path = st.text_input("CSV path",
                              value=DEFAULT_CSV[name], key="csv_path")

    st.markdown("**Computation parameters**")
    max_length = st.slider("Max geodesic length ℓ_max", 4.0, 12.0, 6.0, 0.5)
    mn_max     = st.slider("Truncation m+n ≤ MN_MAX", 2, 10, 4)
    max_prims  = st.slider("Max primitive geodesics", 200, 5000, 1000, 100)
    r_points   = st.slider("Scan points", 100, 1000, 300, 50)
    r_max      = st.slider("r_max", 3.0, 10.0, 6.0, 0.5)
    n_eigs     = st.slider("Eigenvalues to find", 1, 8, 5)

    st.markdown("**Checkpointing**")
    ckpt_dir = st.text_input("Checkpoint directory",
                              value=DEFAULT_CHECKPOINT_DIR)
    use_ckpt  = st.checkbox("Resume from checkpoint", value=True)
    save_ckpt = st.checkbox("Save checkpoint", value=True)

    st.markdown("---")
    run_btn = st.button("▶  Run computation", type="primary",
                        use_container_width=True)
    st.markdown("---")
    st.caption(f"vol({name}) = {VOLUMES[name]:.4f}  |  "
               f"H₁ = ℤ/5  |  ref λ₁ ≈ {REF_LAMBDA[name]}")


# ── Main panel ─────────────────────────────────────────────────────
st.title("Selberg Zeta Eigenvalue Computation")
st.markdown(
    "Estimates Laplace eigenvalues $\\lambda_k = 1 + r_k^2$ for compact "
    "hyperbolic 3-manifolds via the truncated Selberg zeta product."
)

col_info, col_status = st.columns([2, 1])
with col_info:
    st.markdown(f"""
<div class="metric-card">
<b>Manifold:</b> <span class="mono">{name}</span> &nbsp;|&nbsp;
<b>vol:</b> {VOLUMES[name]:.4f} &nbsp;|&nbsp;
<b>H₁:</b> ℤ/5 &nbsp;|&nbsp;
<b>ℓ_max:</b> {max_length} &nbsp;|&nbsp;
<b>MN_MAX:</b> {mn_max} &nbsp;|&nbsp;
<b>r ∈</b> [0.5, {r_max}]
</div>
""", unsafe_allow_html=True)

params = {"name": name, "max_length": max_length, "mn_max": mn_max,
          "max_prims": max_prims, "r_points": r_points, "r_max": r_max}
ckpt_file = checkpoint_path(ckpt_dir, name, params)

with col_status:
    if ckpt_file.exists():
        mtime = datetime.fromtimestamp(ckpt_file.stat().st_mtime)
        st.markdown(
            f'<span class="checkpoint-badge">✓ checkpoint {mtime.strftime("%H:%M %d/%m")}</span>',
            unsafe_allow_html=True)
    else:
        st.markdown('<span style="color:#6b7280;font-size:0.8em">no checkpoint</span>',
                    unsafe_allow_html=True)

# ── Progress area ──────────────────────────────────────────────────
steps_area   = st.empty()
progress_bar = st.empty()
log_area     = st.empty()
results_area = st.empty()

if not run_btn:
    st.markdown("""
<div class="warning-box">
Configure parameters in the sidebar and click <b>▶ Run computation</b>.
<br><br>
<b>Recommended starting point:</b> ℓ_max=6, MN_MAX=4, 1000 primitives → ~60s
<br>
<b>Better accuracy:</b> ℓ_max=10, MN_MAX=6, 3000 primitives → ~15 min
<br>
<b>High accuracy:</b> ℓ_max=12, MN_MAX=8, 5000 primitives → ~1-2 hrs
</div>
""", unsafe_allow_html=True)
    st.stop()

# ── Run ────────────────────────────────────────────────────────────
log_lines = []
def append_log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    log_lines.append(f"[{ts}] {msg}")
    log_area.markdown(
        f'<div class="log-box">' +
        "\n".join(log_lines[-40:]) +
        "</div>", unsafe_allow_html=True)

def show_steps(current):
    steps = [
        ("Load geodesics",      1),
        ("Extract primitives",  2),
        ("Scan log|Z|",         3),
        ("Find minima",         4),
        ("Analyse results",     5),
    ]
    html = ""
    for label, idx in steps:
        if idx < current:
            html += f'<div class="progress-step step-done">✓ {label}</div>'
        elif idx == current:
            html += f'<div class="progress-step step-run">▶ {label}</div>'
        else:
            html += f'<div class="progress-step step-wait">· {label}</div>'
    steps_area.markdown(html, unsafe_allow_html=True)

# ── Check checkpoint ───────────────────────────────────────────────
ckpt_data = None
if use_ckpt:
    ckpt_data = load_checkpoint(ckpt_file)
    if ckpt_data:
        append_log(f"Loaded checkpoint: {ckpt_file}")
        append_log(f"  Saved at: {ckpt_data.get('timestamp', '?')}")
        append_log(f"  Scan points: {len(ckpt_data.get('r_vals', []))}")

# ── Step 1: Load ───────────────────────────────────────────────────
show_steps(1)
append_log(f"Loading geodesics from {csv_path}...")
t0 = time.time()
try:
    ell_all, phi_all = load_geodesics(csv_path, max_length)
    append_log(f"  Loaded {len(ell_all):,} geodesics (ℓ ≤ {max_length})")
except Exception as e:
    st.error(f"Failed to load CSV: {e}")
    st.stop()

# ── Step 2: Primitives ─────────────────────────────────────────────
show_steps(2)
append_log("Extracting primitive geodesics...")
ell_p, phi_p = extract_primitives(ell_all, phi_all, max_prims)
append_log(f"  {len(ell_p)} primitives (max {max_prims})")
append_log(f"  Systole: {ell_p.min():.6f}")

# ── Step 3: Scan (with checkpoint resume) ─────────────────────────
show_steps(3)
r_vals = np.linspace(0.5, r_max, r_points)

if ckpt_data and "log_amp" in ckpt_data:
    log_amp = np.array(ckpt_data["log_amp"])
    r_vals_saved = np.array(ckpt_data["r_vals"])
    if len(r_vals_saved) == len(r_vals) and np.allclose(r_vals_saved, r_vals, atol=1e-6):
        append_log("  Using cached scan from checkpoint — skipping computation")
        r_vals = r_vals_saved
    else:
        append_log("  Checkpoint scan resolution differs — recomputing")
        ckpt_data = None
        log_amp = None
else:
    log_amp = None

if log_amp is None:
    pbar = progress_bar.progress(0, text="Scanning r values...")
    scan_start = time.time()

    def prog_cb(frac, elapsed, eta):
        pbar.progress(min(frac, 1.0),
                      text=f"Scanning… {frac*100:.0f}%  "
                           f"elapsed {elapsed:.0f}s  ETA {eta:.0f}s")

    def log_cb(msg):
        append_log(msg)

    log_amp = run_scan(ell_p, phi_p, mn_max, r_vals, prog_cb, log_cb)
    pbar.progress(1.0, text="Scan complete.")
    append_log(f"Scan complete in {time.time()-scan_start:.1f}s")

    # Save checkpoint
    if save_ckpt:
        ckpt_data_new = {
            "timestamp":  datetime.now().isoformat(),
            "params":     params,
            "r_vals":     r_vals.tolist(),
            "log_amp":    log_amp.tolist(),
            "n_prims":    len(ell_p),
        }
        save_checkpoint(ckpt_file, ckpt_data_new)
        append_log(f"Checkpoint saved: {ckpt_file}")

# ── Step 4: Find minima ────────────────────────────────────────────
show_steps(4)
append_log("Finding eigenvalues from minima of log|Z|...")
minima = find_minima(r_vals, log_amp, ell_p, phi_p, mn_max, n_eigs)
append_log(f"  Found {len(minima)} candidate eigenvalues")
for i, m in enumerate(minima):
    append_log(f"  λ_{i+1} = {m['lambda']:.5f}  (r={m['r']:.5f}, "
               f"log|Z|_min={m['logZ_min']:.3f})")

# ── Step 5: Analyse ────────────────────────────────────────────────
show_steps(5)
append_log("Checking φ-lattice alignment of eigenvalues...")
for m in minima:
    lc = phi_lattice_check(m["lambda"])
    m["lattice"] = lc
    append_log(f"  λ={m['lambda']:.4f}  q*={lc['q']}  "
               f"δ={lc['residual']:.4f}  |Δ|={lc['pct_err']:.1f}%")

append_log(f"Total time: {time.time()-t0:.1f}s")
show_steps(6)  # all done

# ── Results display ────────────────────────────────────────────────
progress_bar.empty()

st.markdown("---")
st.subheader(f"Results — {name}")

# Eigenvalue cards
if minima:
    cols = st.columns(min(len(minima), 5))
    ref = REF_LAMBDA[name]
    for i, (col, m) in enumerate(zip(cols, minima)):
        with col:
            delta = abs(m["lambda"] - ref) / ref * 100 if i == 0 else None
            delta_str = f"<div class='sub'>Δ from ref: {delta:.1f}%</div>" \
                        if delta is not None else ""
            lc = m.get("lattice", {})
            st.markdown(f"""
<div class="eigenvalue-card">
  <div class="sub">λ_{i+1}</div>
  <div class="lam">{m["lambda"]:.4f}</div>
  <div class="sub">r = {m["r"]:.4f}</div>
  {delta_str}
  <div class="sub" style="margin-top:6px;color:#81c784">
    q* = {lc.get("q","?")} &nbsp; δ = {lc.get("residual",0):.3f}
  </div>
</div>
""", unsafe_allow_html=True)

# Reference comparison
st.markdown("---")
c1, c2, c3 = st.columns(3)
with c1:
    st.metric("Reference λ₁ (twist paper)", f"{REF_LAMBDA[name]:.3f}")
with c2:
    if minima:
        st.metric("Computed λ₁", f"{minima[0]['lambda']:.4f}",
                  delta=f"{(minima[0]['lambda']-REF_LAMBDA[name]):+.4f}")
with c3:
    st.metric("log|Z| at minimum", f"{minima[0]['logZ_min']:.3f}" if minima else "—")

# Plot
st.markdown("---")
st.subheader("Selberg zeta scan")
buf = make_plot(name, r_vals, log_amp, minima, REF_LAMBDA[name])
st.image(buf, use_container_width=True)

# φ-lattice table
st.markdown("---")
st.subheader("φ-lattice alignment check")
st.markdown(
    "Are the Laplace eigenvalues near the same lattice "
    r"$\lambda \approx \varphi^{q/4}$ as the particle masses?"
)
if minima:
    rows = []
    for i, m in enumerate(minima):
        lc = m.get("lattice", {})
        rows.append({
            "k": i + 1,
            "λ_k": f"{m['lambda']:.5f}",
            "r_k": f"{m['r']:.5f}",
            "q*": lc.get("q", "?"),
            "φ^(q*/4)": f"{lc.get('predicted', 0):.5f}",
            "Residual δ": f"{lc.get('residual', 0):.4f}",
            "|Δ|/λ (%)": f"{lc.get('pct_err', 0):.2f}%",
        })
    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)

# Download
st.markdown("---")
result_json = {
    "manifold": name, "timestamp": datetime.now().isoformat(),
    "params": params,
    "eigenvalues": [{"k": i+1, **m} for i, m in enumerate(minima)],
}
st.download_button("⬇ Download results JSON",
                   json.dumps(result_json, indent=2),
                   f"selberg_{name}_results.json", "application/json")

# Caveats
st.markdown("---")
st.markdown("""
<div class="warning-box">
<b>Important caveats:</b><br>
1. The truncated Selberg product does not converge on the critical line Re(s)=1.
   These are estimates only; accuracy improves with larger ℓ_max and MN_MAX but
   does not converge monotonically.<br>
2. The minimum of log|Z| approximates a true zero of Z(s); it is not a rigorous
   eigenvalue. A true zero requires log|Z| → −∞.<br>
3. Reference values λ₁(m003)≈2.48 and λ₁(m006)≈2.82 from DQ14199 are more
   reliable and should be cited in papers.<br>
4. The φ-lattice check on eigenvalues is exploratory and speculative.
</div>
""", unsafe_allow_html=True)
