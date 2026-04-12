from snappy import *
import numpy as np
from scipy.linalg import logm, qr

def matrix_to_axis_vector(matrix):
    mat = np.array(matrix, dtype=complex)
    mat = mat / np.sqrt(np.linalg.det(mat))
    L = logm(mat)
    x = float(np.real(L[0,1]+L[1,0]))/2
    y = float(np.imag(L[1,0]-L[0,1]))/2
    z = float(np.real(L[0,0]-L[1,1]))/2
    v = np.array([x,y,z])
    n = np.linalg.norm(v)
    return v/n if n>1e-10 else np.array([1.,0.,0.])

def vector_angle(v1,v2):
    return np.arccos(np.clip(np.abs(np.dot(v1,v2)),-1.,1.))

def fitness(rho, words, sigma, target):
    axes = [matrix_to_axis_vector(rho(w)) for w in words]
    mo = np.zeros((3,3),dtype=complex)
    for r in range(3):
        for c in range(3):
            mo[r,c] = np.exp(-vector_angle(axes[r],axes[c])**2/(2*sigma**2))
    U,_ = qr(mo)
    return np.linalg.norm(np.abs(U)-target,'fro')

CKM = np.array([[0.97427,0.22536,0.00355],
                [0.22522,0.97339,0.04108],
                [0.00886,0.04050,0.99914]])

M = OrientableClosedCensus[43]
rho = M.polished_holonomy()
words = ['aaB','AbA','AAb']

# Get first 5 geodesics
ls = M.length_spectrum(3.0)
geos = list(ls)

PHI = (1+5**0.5)/2
log_phi = np.log(PHI)

print("=== m006 length spectrum vs sigma candidates ===")
print()
print(f"{'idx':>3}  {'Re(l)':>10}  {'fitness':>10}  {'ratio/opt':>10}")
print("-"*45)
for i, g in enumerate(geos[:8]):
    l = float(g['length'].real())
    f = fitness(rho, words, l, CKM)
    print(f"{i:3d}  {l:10.6f}  {f:10.6f}  {l/0.488:10.4f}")

print()
print(f"log(phi) = {log_phi:.6f}  fitness = {fitness(rho,words,log_phi,CKM):.6f}")
print(f"sigma_opt= 0.488000  fitness = {fitness(rho,words,0.488,CKM):.6f}")
print()

# The second geodesic
l2 = float(geos[1]['length'].real())
f_l2 = fitness(rho, words, l2, CKM)
print(f"=== Second geodesic as sigma ===")
print(f"l2 = {l2:.6f}")
print(f"fitness(l2) = {f_l2:.6f}")
print(f"fitness(log_phi) = {fitness(rho,words,log_phi,CKM):.6f}")
print(f"fitness(0.488) = {fitness(rho,words,0.488,CKM):.6f}")
print()
print(f"l2 / log_phi = {l2/log_phi:.6f}")
print(f"l2 / sigma_opt = {l2/0.488:.6f}")
print()

# Physical interpretation
print("=== Interpretation ===")
print(f"The second geodesic of m006 has length l2 = {l2:.6f} rad")
print(f"sigma_opt = {0.488:.6f} rad")
print(f"Difference = {abs(l2-0.488):.6f} rad = {np.degrees(abs(l2-0.488)):.3f} deg")
print(f"Ratio l2/sigma_opt = {l2/0.488:.6f}")
print()
print("Candidate derivation:")
print("sigma = l2(m006) = length of second-shortest geodesic of CKM manifold")
print(f"This would make sigma a spectral property of m006, not a free parameter.")
print()
print(f"Compare: log(phi) = {log_phi:.6f}, l2 = {l2:.6f}")
print(f"l2 - log(phi) = {l2-log_phi:.6f} rad = {np.degrees(l2-log_phi):.3f} deg")
print(f"Is l2 ~ log(phi)? ratio = {l2/log_phi:.6f}")
