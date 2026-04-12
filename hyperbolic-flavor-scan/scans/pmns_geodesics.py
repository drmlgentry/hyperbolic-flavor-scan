import snappy
import numpy as np
from scipy.linalg import logm, qr
from itertools import combinations

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
    for i in range(3):
        for j in range(3):
            mo[i,j] = np.exp(-vector_angle(axes[i],axes[j])**2/(2*sigma**2))
    U,_ = qr(mo)
    return np.linalg.norm(np.abs(U)-target,'fro')

PMNS = np.array([[0.822,0.547,0.152],
                 [0.430,0.611,0.664],
                 [0.371,0.572,0.732]])

M_p = snappy.OrientableClosedCensus[1]
M_c = snappy.OrientableClosedCensus[43]
rho_p = M_p.polished_holonomy()
rho_c = M_c.polished_holonomy()

# Geodesics
spec_p = list(M_p.length_spectrum(2.0))
spec_c = list(M_c.length_spectrum(2.0))

l_p = [float(g['length'].real()) for g in spec_p[:5]]
l_c = [float(g['length'].real()) for g in spec_c[:5]]

print("=== m003 (PMNS) geodesics ===")
for i,l in enumerate(l_p):
    print(f"  l{i+1} = {l:.6f}")

print()
print("=== m006 (CKM) geodesics ===")
for i,l in enumerate(l_c):
    print(f"  l{i+1} = {l:.6f}")

print()
print("=== CKM: sigma_opt=0.488 vs l2=0.491 ===")
print(f"  l2(m006) = {l_c[1]:.6f}  ratio to sigma_opt = {l_c[1]/0.488:.4f}")

print()
print("=== PMNS: sigma_opt=0.500 vs geodesics ===")
sigma_pmns = 0.500
for i,l in enumerate(l_p[:5]):
    print(f"  l{i+1}(m003) = {l:.6f}  ratio to sigma_opt = {l/sigma_pmns:.4f}")

print()
# Key question: is sigma_PMNS / sigma_CKM related to a manifold ratio?
print("=== Ratio of optimal sigmas ===")
print(f"sigma_PMNS / sigma_CKM = {sigma_pmns/0.488:.6f}")
print(f"l1(m003) / l1(m006)    = {l_p[0]/l_c[0]:.6f}")
print(f"l2(m003) / l2(m006)    = {l_p[1]/l_c[1]:.6f}")
print(f"vol(m003) / vol(m006)  = {float(M_p.volume())/float(M_c.volume()):.6f}")
print(f"sqrt(vol ratio)        = {(float(M_p.volume())/float(M_c.volume()))**0.5:.6f}")
print()

# Best PMNS words
words_pmns = ['ab','bA','bba']
words_ckm  = ['aaB','AbA','AAb']

print("=== PMNS axis angles ===")
axes_p = [matrix_to_axis_vector(rho_p(w)) for w in words_pmns]
for i in range(3):
    for j in range(i+1,3):
        a = vector_angle(axes_p[i], axes_p[j])
        print(f"  {words_pmns[i]}-{words_pmns[j]}: {np.degrees(a):.3f} deg")

print()
print("=== CKM axis angles ===")
axes_c = [matrix_to_axis_vector(rho_c(w)) for w in words_ckm]
for i in range(3):
    for j in range(i+1,3):
        a = vector_angle(axes_c[i], axes_c[j])
        print(f"  {words_ckm[i]}-{words_ckm[j]}: {np.degrees(a):.3f} deg")

print()
print("=== The l2 pattern ===")
print(f"CKM: l2 = {l_c[1]:.4f}, sigma_opt = 0.488, ratio = {l_c[1]/0.488:.4f}")
print(f"PMNS: l2 = {l_p[1]:.4f}, sigma_opt = 0.500, ratio = {l_p[1]/0.500:.4f}")
print()
print(f"PMNS fitness at l2(m003) = {fitness(rho_p, words_pmns, l_p[1], PMNS):.5f}")
print(f"PMNS fitness at sigma=0.500 = {fitness(rho_p, words_pmns, 0.500, PMNS):.5f}")
print()
print("Does l2 pattern hold for PMNS?")
print(f"  l2(m003)/sigma_PMNS = {l_p[1]/sigma_pmns:.4f}")
print(f"  l2(m006)/sigma_CKM  = {l_c[1]/0.488:.4f}")
