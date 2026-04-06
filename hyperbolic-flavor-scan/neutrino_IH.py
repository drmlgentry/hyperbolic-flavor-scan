import numpy as np

PHI = (1+5**0.5)/2
m_e = 0.51099895e6  # eV

# PDG inverted ordering
Dm21_sq_PDG = 7.53e-5
Dm32_sq_PDG = 2.51e-3  # |Dm32| for IH ~ same magnitude
sigma_Dm21  = 0.18e-5
sigma_Dm31  = 0.03e-3

def chi2_IH(Dm21, Dm32):
    return ((Dm21-Dm21_sq_PDG)/sigma_Dm21)**2 + \
           ((abs(Dm32)-Dm32_sq_PDG)/sigma_Dm31)**2

print("=== INVERTED ORDERING SCAN ===")
print("IH: m3 < m1 < m2, Dm32^2 < 0")
print()

solutions = []
for q3 in range(-300, -50):
    m3 = m_e * PHI**(q3/4)*1e3
    if not (0.1 < m3 < 30): continue
    for q1 in range(q3+1, q3+80):
        m1 = m_e * PHI**(q1/4)*1e3
        if m1 <= m3: continue
        Dm31 = (m3*1e-3)**2 - (m1*1e-3)**2  # negative
        for q2 in range(q1+1, q1+40):
            m2 = m_e * PHI**(q2/4)*1e3
            if m2 <= m1: continue
            Dm21 = (m2*1e-3)**2 - (m1*1e-3)**2
            if not (1e-5 < Dm21 < 3e-4): continue
            Dm32 = (m3*1e-3)**2 - (m2*1e-3)**2  # negative
            if not (1e-3 < abs(Dm32) < 5e-3): continue
            c2 = chi2_IH(Dm21, Dm32)
            if c2 < 9.0:
                solutions.append(
                    (c2,q3,q1,q2,m3,m1,m2,Dm21,Dm32))

solutions.sort()
print(f"IH solutions within 1-sigma: {sum(1 for s in solutions if s[0]<1)}")
print(f"IH solutions within 3-sigma: {len(solutions)}")
print()
if solutions:
    print(f"{'chi2':>8} {'q3':>5} {'q1':>5} {'q2':>5} "
          f"{'m3':>8} {'m1':>8} {'m2':>8} {'sum':>8}")
    print("-"*65)
    for s in solutions[:5]:
        c2,q3,q1,q2,m3,m1,m2,d21,d32 = s
        print(f"{c2:>8.4f} {q3:>5} {q1:>5} {q2:>5} "
              f"{m3:>8.3f} {m1:>8.3f} {m2:>8.3f} "
              f"{m1+m2+m3:>8.1f}")
else:
    print("No IH solutions found within 3-sigma.")
    print("This would mean the lattice predicts normal ordering.")
