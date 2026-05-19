from snappy import *
import numpy as np
from scipy.linalg import logm, qr

PHI = (1+5**0.5)/2
log_phi = np.log(PHI)

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

M_ckm = OrientableClosedCensus[43]
rho_ckm = M_ckm.polished_holonomy()
words = ['aaB','AbA','AAb']

print("=== Candidate geometric derivations of sigma ===")
print()

# Get manifold invariants
ls = M_ckm.length_spectrum(2.0)
geos = list(ls)
l1 = float(geos[0]['length'].real())
l2 = float(geos[1]['length'].real()) if len(geos)>1 else None
print(f"m006 shortest geodesic: l1 = {l1:.6f}")
print(f"m006 second geodesic:   l2 = {l2:.6f}" if l2 else "")

# CS invariant
cs = float(Manifold('m006').chern_simons())
print(f"CS(m006 cusped) = {cs:.6f}")
print(f"|CS| = {abs(cs):.6f}")
print()

# Candidate sigma values and their fitness
candidates = {
    'log(phi)':          log_phi,
    'l1/pi':             l1/np.pi,
    'l1/2':              l1/2,
    'l1*log(phi)':       l1*log_phi,
    'sqrt(l1*log(phi))': np.sqrt(l1*log_phi),
    'log(phi)/l1':       log_phi/l1,
    '1/phi':             1/PHI,
    '1/phi^2':           1/PHI**2,
    'phi-1':             PHI-1,
    'log(2)':            np.log(2),
    '|CS|':              abs(cs),
    'pi/l1^2':           np.pi/l1**2,
    '2*log(phi)':        2*log_phi,
    'sigma_opt':         0.488,
}

print(f"{'Candidate':30s}  {'value':8s}  {'fitness':10s}  {'ratio/opt':10s}")
print("-"*65)
for name, val in sorted(candidates.items(), key=lambda x: 
                         fitness(rho_ckm, words, x[1], CKM)):
    f = fitness(rho_ckm, words, val, CKM)
    print(f"{name:30s}  {val:.5f}  {f:.6f}  {val/0.488:.4f}")

print()

# Fine scan: which formula matches sigma_opt=0.488 best?
print("=== Checking l1-based formulas ===")
print(f"l1 = {l1:.6f}")
print(f"sigma_opt = 0.488000")
print(f"sigma_opt / l1 = {0.488/l1:.6f}")
print(f"sigma_opt * l1 = {0.488*l1:.6f}")
print(f"sigma_opt^2 / l1 = {0.488**2/l1:.6f}")
print(f"l1 / sigma_opt = {l1/0.488:.6f}")
print()

# The pairwise angles between CKM axes
axes = [matrix_to_axis_vector(rho_ckm(w)) for w in words]
angles = []
for i in range(3):
    for j in range(i+1,3):
        a = vector_angle(axes[i], axes[j])
        angles.append(a)
        print(f"Angle {words[i]}-{words[j]}: {np.degrees(a):.4f} deg = {a:.6f} rad")

print()
min_angle = min(angles)
max_angle = max(angles)
mean_angle = np.mean(angles)
print(f"Min angle:  {np.degrees(min_angle):.4f} deg = {min_angle:.6f} rad")
print(f"Max angle:  {np.degrees(max_angle):.4f} deg = {max_angle:.6f} rad")
print(f"Mean angle: {np.degrees(mean_angle):.4f} deg = {mean_angle:.6f} rad")
print()
print(f"sigma_opt / min_angle  = {0.488/min_angle:.4f}")
print(f"sigma_opt / mean_angle = {0.488/mean_angle:.4f}")
print(f"sigma_opt / max_angle  = {0.488/max_angle:.4f}")
print()
print(f"=== Key ratios ===")
print(f"sigma_opt = {0.488:.5f}")
print(f"log(phi)  = {log_phi:.5f}  ratio = {0.488/log_phi:.5f}")
print(f"l1/pi     = {l1/np.pi:.5f}  ratio = {0.488/(l1/np.pi):.5f}")
print(f"mean_ang  = {mean_angle:.5f}  ratio = {0.488/mean_angle:.5f}")
