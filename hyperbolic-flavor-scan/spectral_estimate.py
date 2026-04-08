import snappy
import numpy as np

PHI = (1+5**0.5)/2

print("=== Spectral gap estimates via length spectrum ===")
print("Using Selberg trace formula lower bound:")
print("lambda_1 >= 1/4 + (pi/l_1)^2")
print()

for idx, name in [(0,'Weeks'),(1,'PMNS'),(43,'CKM')]:
    M = snappy.OrientableClosedCensus[idx]
    ls = M.length_spectrum(3.0)
    geos = list(ls)
    
    print(f"{name} (vol={float(M.volume()):.4f}):")
    
    # Spectral estimates from first few geodesics
    for i, g in enumerate(geos[:5]):
        l = float(g['length'].real())
        # Huber/Selberg: each geodesic of length l contributes
        # to spectrum near (1/4 + (2*pi*n/l)^2) for n=1,2,...
        for n in range(1, 4):
            lam = 0.25 + (2*np.pi*n/l)**2
            # Is this near a phi-lattice point?
            log_lam = np.log(lam)/np.log(PHI)
            nearest_q = round(log_lam/0.25)*0.25
            residual = abs(log_lam - nearest_q)
            marker = " <-- phi-lattice" if residual < 0.05 else ""
            if i < 2:  # show first 2 geodesics
                print(f"  geo{i+1} l={l:.4f} n={n}: "
                      f"lambda~{lam:.4f} "
                      f"log_phi={log_lam:.3f} "
                      f"resid={residual:.3f}{marker}")
    print()
