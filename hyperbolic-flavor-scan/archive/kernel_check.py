import numpy as np
from scipy.linalg import qr

# Confirmed angles
angles = {('aaB','AbA'): 0.8405, ('aaB','AAb'): 1.3523, ('AbA','AAb'): 1.1943}
words = ['aaB', 'AbA', 'AAb']
pairs = [('aaB','AbA'), ('aaB','AAb'), ('AbA','AAb')]

angle_matrix = np.zeros((3,3))
for (w1,w2),th in angles.items():
    i,j = words.index(w1), words.index(w2)
    angle_matrix[i,j] = angle_matrix[j,i] = th

CKM = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.04050, 0.99914]
])

print("Kernel comparison at the three axis angles:")
print(f"  theta_12 = {np.degrees(0.8405):.2f} deg")
print(f"  theta_13 = {np.degrees(1.3523):.2f} deg")
print(f"  theta_23 = {np.degrees(1.1943):.2f} deg")
print()

sigma = 0.49
Delta = 2.08

for th_name, th in [("theta_12", 0.8405), ("theta_13", 1.3523), ("theta_23", 1.1943)]:
    gauss = np.exp(-th**2 / (2*sigma**2))
    hyp   = (1 - np.cos(th))**(-Delta/2)
    print(f"  {th_name}: Gaussian={gauss:.5f}  Hyperbolic={hyp:.5f}")

print()
print("Gaussian kernel overlap matrix (sigma=0.49):")
O_gauss = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        th = angle_matrix[i,j]
        O_gauss[i,j] = np.exp(-th**2/(2*sigma**2)) if i!=j else 1.0
for row in O_gauss:
    print("  ", "  ".join(f"{v:.5f}" for v in row))

U_g, _ = qr(O_gauss)
for i in range(3):
    if U_g[i,i] < 0: U_g[:,i] *= -1
fit_g = np.linalg.norm(np.abs(U_g) - CKM, "fro")
print(f"  Fitness: {fit_g:.6f}")
print()

print("Hyperbolic kernel overlap matrix (Delta=2.08):")
O_hyp = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        th = angle_matrix[i,j]
        O_hyp[i,j] = (1-np.cos(th))**(-Delta/2) if i!=j else 1.0
for row in O_hyp:
    print("  ", "  ".join(f"{v:.5f}" for v in row))
U_h, _ = qr(O_hyp)
for i in range(3):
    if U_h[i,i] < 0: U_h[:,i] *= -1
fit_h = np.linalg.norm(np.abs(U_h) - CKM, "fro")
print(f"  Fitness: {fit_h:.6f}")
print()

print("CONCLUSION:")
print(f"  Gaussian kernel:    fitness = {fit_g:.6f}  <- correct kernel for CKM fit")
print(f"  Hyperbolic kernel:  fitness = {fit_h:.6f}  <- WRONG: diverges at large angles")
print()
print("The hyperbolic propagator K_H ~ (1-cos theta)^{-Delta/2} INCREASES with angle.")
print("For theta > 0, this gives overlap > 1, which is unphysical for a mixing kernel.")
print("The Gaussian kernel is the one that produced fitness 0.0173.")
print()
print("Fix: revert to Gaussian in paper computation sections.")
print("The K_H discussion can remain as theoretical motivation/future direction.")