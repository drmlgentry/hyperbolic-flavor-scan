"""
baryon_lattice_test.py
Fit baryons to golden ratio lattice and run log-uniform null test.
"""

import numpy as np

phi = (1 + np.sqrt(5)) / 2
m_e = 0.000510999  # GeV
denom = 4

def mass_from_q(q):
    return m_e * phi**(q / denom)

def log_phi_ratio(m):
    return np.log(m / m_e) / np.log(phi)

def best_q(m):
    return int(round(log_phi_ratio(m) * denom))

def rms_residual(masses):
    res = []
    for m in masses:
        q = best_q(m)
        m_pred = mass_from_q(q)
        delta = log_phi_ratio(m) - q/denom
        res.append(delta**2)
    return np.sqrt(np.mean(res))

# Baryon data
baryons = [
    ("p", 0.93827208816),
    ("n", 0.93956542052),
    ("Λ", 1.115683),
    ("Σ⁺", 1.18937),
    ("Σ⁰", 1.192642),
    ("Σ⁻", 1.197449),
    ("Ξ⁰", 1.31486),
    ("Ξ⁻", 1.32171),
    ("Ω⁻", 1.67245),
    ("Δ⁺⁺", 1.232),
]

masses = [m for _, m in baryons]
names = [n for n, _ in baryons]

# Fit and print best q values
print("Baryon fit to golden ratio lattice (m = m_e * φ^(q/4))")
print("=" * 70)
print(f"{'Name':>6}  {'m_obs (GeV)':>12}  {'q':>4}  {'m_pred (GeV)':>12}  {'δ (log_φ)':>10}  {'|Δ|/m (%)':>10}")
print("-" * 70)

results = []
for name, m_obs in baryons:
    q = best_q(m_obs)
    m_pred = mass_from_q(q)
    delta = log_phi_ratio(m_obs) - q/denom
    rel_err = abs(m_pred - m_obs) / m_obs * 100
    results.append((name, m_obs, q, m_pred, delta, rel_err))
    print(f"{name:>6}  {m_obs:12.6f}  {q:4d}  {m_pred:12.6f}  {delta:10.4f}  {rel_err:10.2f}")

obs_rms = rms_residual(masses)
print("-" * 70)
print(f"RMS residual (log_φ units): {obs_rms:.5f}")
print(f"Maximum possible RMS for d={denom}: {0.5/denom:.4f}")
print(f"Fraction of max: {obs_rms / (0.5/denom):.2%}")

# Log-uniform null test
n_trials = 10000
log_min = np.log(min(masses))
log_max = np.log(max(masses))
null_rms = []
for _ in range(n_trials):
    rand_log = np.random.uniform(log_min, log_max, len(masses))
    rand_masses = np.exp(rand_log)
    null_rms.append(rms_residual(rand_masses))

null_rms = np.array(null_rms)
p_value = np.mean(null_rms <= obs_rms)
print(f"\nLog-uniform null test (n={n_trials}): p = {p_value:.4f}")
if p_value < 0.05:
    print("→ Significant clustering at 5% level")
else:
    print("→ Not significant at 5% level")