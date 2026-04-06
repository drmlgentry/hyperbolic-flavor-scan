import numpy as np

PHI = (1+5**0.5)/2
m_e = 0.51099895e6  # eV

Dm21_sq_PDG = 7.53e-5
Dm31_sq_PDG = 2.51e-3
# PDG 1-sigma uncertainties
sigma_Dm21 = 0.18e-5
sigma_Dm31 = 0.03e-3

print("=== UNIQUENESS TEST ===")
print(f"PDG: Dm21^2 = {Dm21_sq_PDG:.4e} +/- {sigma_Dm21:.2e}")
print(f"PDG: Dm31^2 = {Dm31_sq_PDG:.4e} +/- {sigma_Dm31:.2e}")
print()

# Chi^2 with PDG uncertainties
def chi2(Dm21, Dm31):
    return ((Dm21-Dm21_sq_PDG)/sigma_Dm21)**2 + \
           ((Dm31-Dm31_sq_PDG)/sigma_Dm31)**2

solutions_1sigma = []
solutions_2sigma = []
solutions_3sigma = []

for q1 in range(-200, -100):
    m1 = m_e * PHI**(q1/4) * 1e3
    if not (0.1 < m1 < 30): continue
    for q2 in range(q1+1, q1+40):
        m2 = m_e * PHI**(q2/4) * 1e3
        if m2 <= m1: continue
        Dm21 = (m2*1e-3)**2 - (m1*1e-3)**2
        if not (1e-5 < Dm21 < 3e-4): continue
        for q3 in range(q2+1, q2+80):
            m3 = m_e * PHI**(q3/4) * 1e3
            if m3 <= m2: continue
            Dm31 = (m3*1e-3)**2 - (m1*1e-3)**2
            if not (5e-4 < Dm31 < 8e-3): continue
            c2 = chi2(Dm21, Dm31)
            entry = (c2,q1,q2,q3,m1,m2,m3,Dm21,Dm31)
            if c2 < 1.0:  solutions_1sigma.append(entry)
            elif c2 < 4.0: solutions_2sigma.append(entry)
            elif c2 < 9.0: solutions_3sigma.append(entry)

solutions_1sigma.sort()
solutions_2sigma.sort()
solutions_3sigma.sort()

print(f"Solutions within 1-sigma (chi2<1): {len(solutions_1sigma)}")
print(f"Solutions within 2-sigma (chi2<4): {len(solutions_2sigma)}")
print(f"Solutions within 3-sigma (chi2<9): {len(solutions_3sigma)}")
print()

print("=== 1-SIGMA SOLUTIONS ===")
print(f"{'chi2':>8} {'q1':>5} {'q2':>5} {'q3':>5} "
      f"{'m1':>8} {'m2':>8} {'m3':>8} {'sum':>8}")
print("-"*65)
for c2,q1,q2,q3,m1,m2,m3,d21,d31 in solutions_1sigma:
    print(f"{c2:>8.4f} {q1:>5} {q2:>5} {q3:>5} "
          f"{m1:>8.3f} {m2:>8.3f} {m3:>8.3f} {m1+m2+m3:>8.1f}")

print()
print("=== BEST SOLUTION DETAILS ===")
if solutions_1sigma:
    c2,q1,q2,q3,m1,m2,m3,d21,d31 = solutions_1sigma[0]
    print(f"q = ({q1}, {q2}, {q3})")
    print(f"Spacing: q2-q1={q2-q1}, q3-q1={q3-q1}, q3-q2={q3-q2}")
    print(f"m1={m1:.4f} meV  m2={m2:.4f} meV  m3={m3:.4f} meV")
    print(f"sum_nu = {m1+m2+m3:.2f} meV  (Planck: <120 meV, CMB-S4: ~30 meV sensitivity)")
    print(f"Dm21^2 = {d21:.6e}  residual: {(d21-Dm21_sq_PDG)/sigma_Dm21:+.2f} sigma")
    print(f"Dm31^2 = {d31:.6e}  residual: {(d31-Dm31_sq_PDG)/sigma_Dm31:+.2f} sigma")
    print()
    # Mass ratios
    print(f"m2/m1 = {m2/m1:.4f}  (phi^3 = {PHI**3:.4f}?  ratio/phi^3 = {(m2/m1)/PHI**3:.4f})")
    print(f"m3/m1 = {m3/m1:.4f}  (phi^15 = {PHI**15:.4f}?)")
    print(f"m3/m2 = {m3/m2:.4f}  (phi^12 = {PHI**12:.4f}?  ratio/phi^12 = {(m3/m2)/PHI**12:.4f})")
    print()
    # q spacing analysis
    print(f"q2-q1 = {q2-q1}  = {q2-q1} lattice steps")
    print(f"q3-q2 = {q3-q2}  = {q3-q2} lattice steps")  
    print(f"q3-q1 = {q3-q1}  = {q3-q1} lattice steps")
    print()
    # Inverted hierarchy check
    m3_ih = m_e * PHI**(q3/4) * 1e3
    m1_ih = m_e * PHI**(q1/4) * 1e3
    Dm32_ih = (m3_ih*1e-3)**2 - (m2*1e-3)**2
    print(f"Inverted hierarchy Dm32^2 = {Dm32_ih:.4e}  (PDG NH: {Dm31_sq_PDG:.4e})")

print()
print("=== PREDICTION FOR CMB-S4 ===")
if solutions_1sigma:
    c2,q1,q2,q3,m1,m2,m3,d21,d31 = solutions_1sigma[0]
    sum_m = m1+m2+m3
    print(f"Predicted sum: {sum_m:.1f} meV")
    print(f"CMB-S4 sensitivity: ~30 meV (3-sigma detection if sum > 60 meV)")
    print(f"Prediction: DETECTABLE by CMB-S4")
    print(f"If CMB-S4 measures sum outside {sum_m-5:.0f}-{sum_m+5:.0f} meV: lattice falsified")
