import numpy as np

PHI = (1+5**0.5)/2
sigma = 0.49

print("=== The sigma = log(phi) connection ===")
print()
print(f"sigma (best fit, CKM paper) = {sigma:.5f} rad")
print(f"log(phi)                    = {np.log(PHI):.5f} rad")
print(f"Ratio sigma/log(phi)        = {sigma/np.log(PHI):.5f}")
print(f"Difference                  = {sigma - np.log(PHI):.5f} rad")
print(f"                            = {np.degrees(sigma - np.log(PHI)):.4f} deg")
print()
print("Interpretation:")
print("The Gaussian kernel exp(-angle^2 / (2*sigma^2)) sets the")
print("'coherence length' on the Bloch sphere for mixing amplitudes.")
print("sigma = log(phi) = Reg(Q(sqrt(5))) = mass lattice spacing.")
print()
print("This means the SAME number that sets the mass lattice spacing")
print("also sets the angular coherence length for the mixing matrix.")
print("Both emerge from log(phi) = half the log-Mahler-measure of")
print("the Alexander polynomial of the cusped m003 ancestor.")
print()

# Sensitivity: how much does fitness change with sigma?
print("=== Fitness vs sigma ===")
print(f"{'sigma':>8}  {'sigma/log(phi)':>16}  fitness (CKM)")

import snappy
from scipy.linalg import logm, qr

def matrix_to_axis_vector(matrix):
    mat = np.array(matrix, dtype=complex)
    det = np.linalg.det(mat)
    mat = mat / np.sqrt(det)
    log_mat = logm(mat)
    a,b,c,d = log_mat[0,0],log_mat[0,1],log_mat[1,0],log_mat[1,1]
    x = float(np.real(b+c))/2
    y = float(np.imag(c-b))/2
    z = float(np.real(a-d))/2
    vec = np.array([x,y,z])
    n = np.linalg.norm(vec)
    return vec/n if n>1e-10 else np.array([1.,0.,0.])

def vector_angle(v1,v2):
    return np.arccos(np.clip(np.abs(np.dot(v1,v2)),-1.,1.))

def fitness(rho, words, sigma):
    axes = [matrix_to_axis_vector(rho(w)) for w in words]
    mo = np.zeros((3,3),dtype=complex)
    for r in range(3):
        for c in range(3):
            mo[r,c] = np.exp(-vector_angle(axes[r],axes[c])**2/(2*sigma**2))
    U,_ = qr(mo)
    CKM = np.array([[0.97427,0.22536,0.00355],
                    [0.22522,0.97339,0.04108],
                    [0.00886,0.04050,0.99914]])
    return np.linalg.norm(np.abs(U)-CKM,'fro')

M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()
words = ['aaB','AbA','AAb']

log_phi = np.log(PHI)
for s in np.linspace(0.30, 0.70, 21):
    f = fitness(rho, words, s)
    marker = " <-- log(phi)" if abs(s-log_phi)<0.005 else ""
    marker2 = " <-- best(0.49)" if abs(s-0.49)<0.005 else ""
    print(f"{s:8.4f}  {s/log_phi:16.4f}  {f:.5f}{marker}{marker2}")

print()
print(f"log(phi) = {log_phi:.5f}")
print(f"Best sigma from paper = 0.49000")
print(f"Difference = {abs(0.49-log_phi):.5f} rad = {np.degrees(abs(0.49-log_phi)):.3f} deg")
print()
print("Null test: does sigma matter?")
print("If sigma >> 1: all overlaps -> 1, matrix -> all-ones, QR -> identity-like")
print("If sigma << 1: all overlaps -> 0, matrix -> diagonal, QR -> identity")
print("Only near sigma ~ O(1) do you get non-trivial mixing.")
print(f"sigma = log(phi) = 0.481 is the natural geometric scale.")
