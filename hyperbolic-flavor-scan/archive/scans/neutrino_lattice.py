import numpy as np

PHI = (1+5**0.5)/2

# PDG 2024 oscillation data
Dm21_sq_PDG = 7.53e-5   # eV^2
Dm31_sq_PDG = 2.51e-3   # eV^2 (normal ordering)

# Mass formula: m = m_e * phi^(q/4) in eV
# m_e = 0.51099895 MeV = 510998.95 eV
m_e = 0.51099895e6  # eV

print("=== Neutrino mass lattice scan ===")
print("Scanning q values near expected neutrino mass range (meV)")
print()

# Neutrino masses should be in range 1-100 meV
# m = m_e * phi^(q/4) ~ 1e-3 eV means phi^(q/4) ~ 2e-9
# q/4 = log(2e-9)/log(phi) ~ -42
# So q ~ -168 for lightest neutrino

best_chi2 = 1e10
best_triple = None

print(f"{'q1':>5} {'q2':>5} {'q3':>5} {'m1_meV':>10} {'m2_meV':>10} "
      f"{'m3_meV':>10} {'Dm21sq':>12} {'Dm31sq':>12} {'chi2':>10}")
print("-"*90)

for q1 in range(-200, -140):
    m1 = m_e * PHI**(q1/4) * 1e3  # meV
    if not (0.5 < m1 < 20): continue
    
    for q2 in range(q1+1, q1+30):
        m2 = m_e * PHI**(q2/4) * 1e3  # meV
        if m2 <= m1: continue
        
        Dm21_sq = (m2*1e-3)**2 - (m1*1e-3)**2  # eV^2
        if not (0.5e-5 < Dm21_sq < 2e-4): continue
        
        for q3 in range(q2+1, q2+60):
            m3 = m_e * PHI**(q3/4) * 1e3  # meV
            if m3 <= m2: continue
            
            Dm31_sq = (m3*1e-3)**2 - (m1*1e-3)**2  # eV^2
            if not (1e-3 < Dm31_sq < 5e-3): continue
            
            # chi^2 against PDG
            chi2 = ((Dm21_sq - Dm21_sq_PDG)/Dm21_sq_PDG)**2 + \
                   ((Dm31_sq - Dm31_sq_PDG)/Dm31_sq_PDG)**2
            
            if chi2 < 0.1:
                print(f"{q1:>5} {q2:>5} {q3:>5} {m1:>10.3f} {m2:>10.3f} "
                      f"{m3:>10.3f} {Dm21_sq:>12.2e} {Dm31_sq:>12.2e} "
                      f"{chi2:>10.6f}")
            
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_triple = (q1,q2,q3,m1,m2,m3,Dm21_sq,Dm31_sq)

print()
print(f"Best fit: chi2={best_chi2:.6f}")
if best_triple:
    q1,q2,q3,m1,m2,m3,d21,d31 = best_triple
    print(f"  q=({q1},{q2},{q3})")
    print(f"  m=({m1:.3f},{m2:.3f},{m3:.3f}) meV")
    print(f"  sum_m = {(m1+m2+m3):.1f} meV")
    print(f"  Dm21^2 = {d21:.4e} eV^2  (PDG: {Dm21_sq_PDG:.4e})")
    print(f"  Dm31^2 = {d31:.4e} eV^2  (PDG: {Dm31_sq_PDG:.4e})")
    print(f"  Residual Dm21: {(d21-Dm21_sq_PDG)/Dm21_sq_PDG*100:+.2f}%")
    print(f"  Residual Dm31: {(d31-Dm31_sq_PDG)/Dm31_sq_PDG*100:+.2f}%")
    print()
    # Is the solution unique?
    print("Uniqueness: are there other solutions with chi2 < 0.01?")
