import snappy, numpy as np
from scipy.linalg import logm
import warnings; warnings.filterwarnings('ignore')

M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()

def get_axis(word):
    mat = np.array(rho(word), dtype=complex)
    mat /= np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    v = np.array([float(np.real(L[0,1]+L[1,0]))/2,
                  float(np.imag(L[1,0]-L[0,1]))/2,
                  float(np.real(L[0,0]-L[1,1]))/2])
    n = np.linalg.norm(v)
    return v/n if n > 1e-10 else v

words = ['aaB', 'AbA', 'AAb']
axes  = [get_axis(w) for w in words]

print("Axis angles:")
for i in range(3):
    for j in range(i+1,3):
        ang = np.degrees(np.arccos(abs(np.dot(axes[i],axes[j]))))
        print(f"  {words[i]}-{words[j]}: {ang:.4f} deg")

print()
PDG = np.array([[0.97435, 0.22500, 0.00369],
                [0.22486, 0.97349, 0.04182],
                [0.00857, 0.04110, 0.99917]])

print("Fitness vs sigma (canonical triple):")
print(f"{'sigma':>8}  {'fitness':>10}  {'V_us_pred':>10}  {'theta_C_pred':>14}  {'theta_C_err%':>13}")
print("-"*60)
best_fitness = 999; best_sigma = 0; best_U = None
for sigma in np.arange(0.30, 0.80, 0.01):
    O = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            ang = np.degrees(np.arccos(abs(np.dot(axes[i],axes[j]))))
            O[i,j] = np.exp(-np.radians(ang)**2 / (2*sigma**2))
    Q, R = np.linalg.qr(O)
    U = np.abs(Q)
    fitness = np.sqrt(np.sum((np.sort(U.flatten()) - np.sort(PDG.flatten()))**2))
    if fitness < best_fitness:
        best_fitness = fitness; best_sigma = sigma; best_U = U.copy()
    if abs(sigma - 0.49) < 0.005:
        Vus = U[0,1]
        tc_pred = np.degrees(np.arcsin(min(Vus,1)))
        tc_pdg  = np.degrees(np.arcsin(0.22500))
        err = abs(tc_pred - tc_pdg)/tc_pdg*100
        print(f"{sigma:>8.2f}  {fitness:>10.6f}  {Vus:>10.5f}  {tc_pred:>14.4f}  {err:>13.4f}%")

print()
print(f"Best sigma:   {best_sigma:.2f}")
print(f"Best fitness: {best_fitness:.6f}")
print()
print("At best sigma, predicted CKM matrix (absolute values):")
print(np.round(best_U, 5))
print()
print("PDG CKM matrix:")
print(PDG)
print()
Vus_best = best_U[0,1]
tc_pred  = np.degrees(np.arcsin(min(Vus_best,1)))
tc_pdg   = np.degrees(np.arcsin(0.22500))
print(f"V_us predicted: {Vus_best:.5f}  (PDG: 0.22500)")
print(f"Cabibbo angle predicted: {tc_pred:.4f} deg")
print(f"Cabibbo angle PDG:       {tc_pdg:.4f} deg")
print(f"Error: {abs(tc_pred-tc_pdg)/tc_pdg*100:.4f}%")
print()
print("Paper claims: 0.4% error at sigma=0.49")
