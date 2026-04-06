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

def scan_base(base):
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
                    solutions.append(
                        (c2,q1,q2,q3,m1,m2,m3,Dm21,Dm31))
    solutions.sort()
    return solutions

# Fine scan including phi exactly
# Use step 0.005 so phi=1.6180 is hit precisely
print("=== BASE SENSITIVITY SCAN (fine, includes phi exactly) ===")
print(f"phi = {PHI:.6f}")
print()

# Build base list: fine grid + phi exactly
bases = sorted(set(
    list(np.arange(1.50, 1.71, 0.005)) + [PHI]
))

print(f"{'base':>9} {'n_1sig':>7} {'n_2sig':>7} {'n_3sig':>7} "
      f"{'best_chi2':>11} {'best_q':>22}")
print("-"*68)

phi_row = None
for b in bases:
    sols = scan_base(b)
    n1 = sum(1 for s in sols if s[0] < 1.0)
    n2 = sum(1 for s in sols if s[0] < 4.0)
    n3 = len(sols)
    best = sols[0] if sols else None
    best_c2 = f"{best[0]:.4f}" if best else "---"
    best_q  = f"({best[1]},{best[2]},{best[3]})" if best else "---"
    marker  = " <-- PHI" if abs(b-PHI) < 1e-6 else ""
    print(f"{b:>9.5f} {n1:>7} {n2:>7} {n3:>7} "
          f"{best_c2:>11} {best_q:>22}{marker}")
    if abs(b-PHI) < 1e-6:
        phi_row = (n1, n2, n3, best)

print()
print("=== VERDICT ===")
if phi_row:
    n1, n2, n3, best = phi_row
    print(f"phi = {PHI:.6f}")
    print(f"  1-sigma solutions: {n1}")
    print(f"  2-sigma solutions: {n2}")
    print(f"  3-sigma solutions: {n3}")
    if best:
        print(f"  Best chi2: {best[0]:.4f}")
        print(f"  Best q: ({best[1]},{best[2]},{best[3]})")
        print(f"  Masses: {best[4]:.3f}, {best[5]:.3f}, {best[6]:.3f} meV")

print()
# Count how many bases have n_1sig >= 1
# to assess how special phi is
print("=== IS PHI SPECIAL? ===")
print("Bases with at least 1 solution within 1-sigma:")
for b in bases:
    sols = scan_base(b)
    n1 = sum(1 for s in sols if s[0] < 1.0)
    if n1 > 0:
        best = sols[0]
        print(f"  b={b:.5f}: {n1} solutions, "
              f"best chi2={best[0]:.4f} at q=({best[1]},{best[2]},{best[3]})")
