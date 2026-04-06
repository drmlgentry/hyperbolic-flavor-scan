import numpy as np

PHI = (1+5**0.5)/2
m_e = 0.51099895e6  # eV

Dm21_sq_PDG = 7.53e-5
Dm31_sq_PDG = 2.51e-3
sigma_Dm21  = 0.18e-5
sigma_Dm31  = 0.03e-3

def chi2(Dm21, Dm31):
    return ((Dm21-Dm21_sq_PDG)/sigma_Dm21)**2 + \
           ((Dm31-Dm31_sq_PDG)/sigma_Dm31)**2

def scan_base(base, label=""):
    solutions = []
    for q1 in range(-300, -50):
        m1 = m_e * base**(q1/4) * 1e3
        if not (0.1 < m1 < 30): continue
        for q2 in range(q1+1, q1+60):
            m2 = m_e * base**(q2/4) * 1e3
            if m2 <= m1: continue
            Dm21 = (m2*1e-3)**2 - (m1*1e-3)**2
            if not (1e-5 < Dm21 < 3e-4): continue
            for q3 in range(q2+1, q2+100):
                m3 = m_e * base**(q3/4) * 1e3
                if m3 <= m2: continue
                Dm31 = (m3*1e-3)**2 - (m1*1e-3)**2
                if not (5e-4 < Dm31 < 8e-3): continue
                c2 = chi2(Dm21, Dm31)
                if c2 < 9.0:
                    solutions.append((c2,q1,q2,q3,m1,m2,m3,Dm21,Dm31))
    solutions.sort()
    return solutions

# ── 1. Correct mass ratio analysis ───────────────────────────────
print("=== CORRECT MASS RATIOS ===")
print("q=(-149,-146,-134): spacings dq = 3, 12, 15")
print()
q1,q2,q3 = -149,-146,-134
m1 = m_e * PHI**(q1/4)*1e3
m2 = m_e * PHI**(q2/4)*1e3
m3 = m_e * PHI**(q3/4)*1e3

dq21 = q2-q1  # = 3
dq32 = q3-q2  # = 12
dq31 = q3-q1  # = 15

print(f"m2/m1 = phi^(dq21/4) = phi^({dq21}/4) = phi^{dq21/4:.4f}")
print(f"  exact: phi^(3/4) = {PHI**(3/4):.6f}")
print(f"  actual: {m2/m1:.6f}  match: {abs(m2/m1 - PHI**(dq21/4))<1e-6}")
print()
print(f"m3/m2 = phi^(dq32/4) = phi^({dq32}/4) = phi^{dq32/4:.4f} = phi^3")
print(f"  exact: phi^3 = {PHI**3:.6f}")
print(f"  actual: {m3/m2:.6f}  match: {abs(m3/m2 - PHI**3)<1e-6}")
print()
print(f"m3/m1 = phi^(dq31/4) = phi^({dq31}/4) = phi^{dq31/4:.4f} = phi^(15/4)")
print(f"  exact: phi^(15/4) = {PHI**(15/4):.6f}")
print(f"  actual: {m3/m1:.6f}  match: {abs(m3/m1 - PHI**(15/4))<1e-6}")
print()
print("KEY EXACT RELATION: m3/m2 = phi^3  (spacing = 12 quarter-steps = 3 full steps)")
print()

# ── 2. Base sensitivity scan ──────────────────────────────────────
print("=== BASE SENSITIVITY SCAN ===")
print(f"{'base':>7} {'n_1sig':>8} {'n_3sig':>8} {'best_chi2':>12} {'best_q':>20}")
print("-"*60)

bases = np.arange(1.50, 1.72, 0.02)
phi_result = None
for b in bases:
    sols = scan_base(b)
    n1 = sum(1 for s in sols if s[0]<1)
    n3 = len(sols)
    best = sols[0] if sols else None
    best_c2 = f"{best[0]:.4f}" if best else "none"
    best_q = f"({best[1]},{best[2]},{best[3]})" if best else "none"
    marker = " <-- phi" if abs(b-PHI)<0.01 else ""
    print(f"{b:>7.4f} {n1:>8} {n3:>8} {best_c2:>12} {best_q:>20}{marker}")
    if abs(b-PHI) < 0.01:
        phi_result = (n1, n3, best)

print()
print("=== VERDICT ===")
if phi_result:
    n1, n3, best = phi_result
    print(f"phi={PHI:.4f}: {n1} solution(s) within 1-sigma")
    print(f"Is phi a local optimum in chi2? -> check neighbors above")

# ── 3. Proper CMB-S4 prediction with uncertainty ─────────────────
print()
print("=== CMB-S4 PREDICTION WITH PROPAGATED UNCERTAINTY ===")

# Propagate Dm^2 uncertainties to individual masses
# m1^2 = m2^2 - Dm21^2,  m2^2 = m3^2 - Dm32^2
# But on lattice: m_i = m_e * phi^(qi/4), so they are FIXED
# The uncertainty comes only from PDG Dm^2 measurements
# The lattice prediction is a point; the question is whether it
# lies within the PDG error ellipse

# Propagate: sum = m1 + m2 + m3
# d(sum)/d(Dm21^2) and d(sum)/d(Dm31^2) via implicit differentiation
# On lattice the q values are fixed integers, so masses are exact
# The chi^2 = 0.58 tells us we are 0.76-sigma from PDG center

c2_best = 0.5803
sigma_sum = 0.76  # sqrt(chi2) in sigma units from PDG center

# More careful: propagate Dm^2 uncertainties to mass sum
# m3^2 = m1^2 + Dm31^2, m2^2 = m1^2 + Dm21^2
# sum = m1 + sqrt(m1^2+Dm21^2) + sqrt(m1^2+Dm31^2)
# On lattice m1 is fixed, vary Dm^2 within 1-sigma

sum_samples = []
for _ in range(100000):
    d21 = np.random.normal(Dm21_sq_PDG, sigma_Dm21)
    d31 = np.random.normal(Dm31_sq_PDG, sigma_Dm31)
    if d21 <= 0 or d31 <= d21: continue
    # Find closest lattice solution
    # Use our best-fit q values, compute what sum would be
    # if Dm^2 were exactly these values (off-lattice masses)
    m1_true = np.sqrt(max(0, m3**2 - d31*1e6))  # meV
    m2_true = np.sqrt(max(0, m1_true**2 + d21*1e6))  # meV  
    m3_true = np.sqrt(max(0, m1_true**2 + d31*1e6))  # meV
    sum_samples.append(m1_true + m2_true + m3_true)

sum_samples = np.array(sum_samples)
print(f"Lattice prediction: sum_nu = {m1+m2+m3:.2f} meV")
print(f"Propagated 1-sigma from PDG Dm^2 uncertainties:")
print(f"  sum = {np.mean(sum_samples):.2f} +/- {np.std(sum_samples):.2f} meV")
print(f"  68% CI: [{np.percentile(sum_samples,16):.1f}, {np.percentile(sum_samples,84):.1f}] meV")
print(f"  95% CI: [{np.percentile(sum_samples,2.5):.1f}, {np.percentile(sum_samples,97.5):.1f}] meV")
print()
print(f"Lattice falsified at 3-sigma if CMB-S4 measures:")
pred = m1+m2+m3
sigma_pred = np.std(sum_samples)
print(f"  sum_nu < {pred - 3*sigma_pred:.1f} meV  or  sum_nu > {pred + 3*sigma_pred:.1f} meV")
print(f"CMB-S4 sensitivity ~30 meV -> will detect sum={pred:.1f} meV at >{pred/30:.1f}-sigma")
print(f"Prediction is TESTABLE and FALSIFIABLE by CMB-S4.")
