"""
lattice_fit_null_test.py
Fit baryons and mesons to golden ratio lattice and run log-uniform null test.
"""

import numpy as np

# Constants
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
    """RMS of residuals (log_phi units) for given list of masses."""
    res = []
    for m in masses:
        q = best_q(m)
        m_pred = mass_from_q(q)
        delta = log_phi_ratio(m) - q/denom
        res.append(delta**2)
    return np.sqrt(np.mean(res))

# ------------------------------------------------------------
# Data: name, mass_GeV (PDG 2024 central values)
# Baryons (already tested)
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
# Mesons (light unflavored and strange, plus charmonium, bottomonium, etc.)
mesons = [
    ("π⁰", 0.1349768),
    ("π⁺", 0.13957039),
    ("π⁻", 0.13957039),
    ("K⁺", 0.493677),
    ("K⁻", 0.493677),
    ("K⁰", 0.497611),
    ("η", 0.547862),
    ("η'", 0.95778),
    ("ρ⁰", 0.77526),
    ("ρ⁺", 0.77526),
    ("ρ⁻", 0.77526),
    ("ω", 0.78265),
    ("φ", 1.019461),
    ("J/ψ", 3.0969),
    ("ψ(2S)", 3.6861),
    ("Υ(1S)", 9.4603),
    ("Υ(2S)", 10.0233),
    ("Υ(3S)", 10.3552),
    ("B⁺", 5.27934),
    ("B⁻", 5.27934),
    ("B⁰", 5.27965),
    ("B_s", 5.36692),
    ("D⁺", 1.86966),
    ("D⁻", 1.86966),
    ("D⁰", 1.86484),
    ("D_s", 1.96834),
    ("η_c", 2.9839),
    ("χ_c0", 3.41475),
    ("χ_c1", 3.51066),
    ("χ_c2", 3.5562),
    ("h_c", 3.52538),
    ("ψ(3770)", 3.7737),
]

# Combine all particles
all_particles = baryons + mesons
masses = [m for _, m in all_particles]
names = [n for n, _ in all_particles]

print("Fitting to golden ratio lattice (m_e * φ^(q/4))")
print(f"Number of particles: {len(masses)}")
print(f"Mass range: {min(masses):.6f} - {max(masses):.6f} GeV")
print()

obs_rms = rms_residual(masses)
print(f"Observed RMS residual: {obs_rms:.5f} (log_φ units)")

# Log-uniform null test
n_trials = 5000
log_min = np.log(min(masses))
log_max = np.log(max(masses))
null_rms = []
for _ in range(n_trials):
    rand_log = np.random.uniform(log_min, log_max, len(masses))
    rand_masses = np.exp(rand_log)
    null_rms.append(rms_residual(rand_masses))

null_rms = np.array(null_rms)
p_value = np.mean(null_rms <= obs_rms)
print(f"Log-uniform null test (n={n_trials}): p = {p_value:.4f}")

if p_value < 0.05:
    print("→ Significant clustering at 5% level")
else:
    print("→ Not significant at 5% level (p > 0.05)")

print(f"Null RMS mean: {null_rms.mean():.5f}  std: {null_rms.std():.5f}")