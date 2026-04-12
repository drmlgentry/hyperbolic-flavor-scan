import snappy
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
    for i in range(3):
        for j in range(3):
            ang = vector_angle(axes[i], axes[j])
            mo[i,j] = np.exp(-ang**2/(2*sigma**2))
    U,_ = qr(mo)
    return np.linalg.norm(np.abs(U)-target,'fro')

PMNS = np.array([[0.822,0.547,0.152],
                 [0.430,0.611,0.664],
                 [0.371,0.572,0.732]])

M = snappy.OrientableClosedCensus[1]
rho = M.polished_holonomy()
words = ['ab','bA','bba']

sigma_vals = np.linspace(0.4,0.6,21)
best_fit = float('inf')
best_sigma = None
for s in sigma_vals:
    f = fitness(rho, words, s, PMNS)
    if f < best_fit:
        best_fit = f
        best_sigma = s
    print(f"{s:.3f} {f:.5f}")
print(f"Best sigma: {best_sigma:.3f} fitness: {best_fit:.5f}")
