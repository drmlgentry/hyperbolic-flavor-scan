import snappy, numpy as np
from scipy.linalg import logm

def matrix_to_axis_vector(matrix):
    mat = np.array(matrix, dtype=complex)
    det = np.linalg.det(mat)
    mat = mat / np.sqrt(det)
    L = logm(mat)
    a, b, c, d = L[0,0], L[0,1], L[1,0], L[1,1]
    x = float(np.real(b + c)) / 2
    y = float(np.imag(c - b)) / 2
    z = float(np.real(a - d)) / 2
    vec = np.array([x, y, z])
    return vec / np.linalg.norm(vec)

M = snappy.OrientableClosedCensus[43]
rho = M.polished_holonomy()
words = ['aaB', 'AbA', 'AAb']
axes = [matrix_to_axis_vector(rho(w)) for w in words]
pairs = [('aaB','AbA',0,1), ('aaB','AAb',0,2), ('AbA','AAb',1,2)]
for w1, w2, i, j in pairs:
    dot = np.clip(np.abs(np.dot(axes[i], axes[j])), -1, 1)
    angle_rad = np.arccos(dot)
    print(f"{w1} <-> {w2}: {angle_rad:.4f} rad = {np.degrees(angle_rad):.1f} deg")
