import snappy
import numpy as np
from scipy.linalg import qr, logm
import itertools

# The winning manifold
M = snappy.OrientableClosedCensus[43]
print("Manifold: " + M.name())
print("Volume: " + str(float(M.volume())))
print("Homology: " + str(M.homology()))

# Get the holonomy
G = M.fundamental_group()
rho = M.polished_holonomy()
print("Generators: " + str(G.generators()))

# The winning words
words = ['aaB', 'aBa', 'AAb']

def matrix_to_axis_vector(matrix):
    mat = np.array(matrix, dtype=complex)
    det = np.linalg.det(mat)
    mat = mat / np.sqrt(det)
    log_mat = logm(mat)
    a, b, c, d = log_mat[0,0], log_mat[0,1], log_mat[1,0], log_mat[1,1]
    x = float(np.real(b + c)) / 2
    y = float(np.imag(c - b)) / 2
    z = float(np.real(a - d)) / 2
    vec = np.array([x, y, z])
    norm = np.linalg.norm(vec)
    return vec / norm if norm > 1e-10 else np.array([1., 0., 0.])

def get_word_homology_image(word_str, num_gens):
    image = np.zeros(num_gens)
    for char in word_str:
        if char.islower():
            idx = ord(char) - ord('a')
            if idx < num_gens: image[idx] += 1
        elif char.isupper():
            idx = ord(char) - ord('A')
            if idx < num_gens: image[idx] -= 1
    return image

# Compute axes
axes = []
for w in words:
    mat = rho(w)
    axes.append(matrix_to_axis_vector(mat))
    print("Word '" + w + "': axis = " + str(np.round(axes[-1], 4)))

# Compute homology images
print("\nHomology images:")
for w in words:
    img = get_word_homology_image(w, 2)
    print("  " + w + " -> " + str(img))

# Now scan over CP phases
print("\n" + "="*60)
print("Scanning CP phases to match Jarlskog")
print("="*60)

V_CKM_TARGET = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.0405, 0.99914]
])
J_TARGET = 3.0e-5
sigma = 0.5

def vector_angle(v1, v2):
    dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
    return np.arccos(np.abs(dot))

best_result = None
best_j_error = float('inf')

# Scan over phase angles
for phi1 in np.linspace(0, 2*np.pi, 20):
    for phi2 in np.linspace(0, 2*np.pi, 20):
        char_vec = np.array([phi1, phi2])
        
        phases = []
        for w in words:
            hom_img = get_word_homology_image(w, 2)
            phase = np.dot(hom_img, char_vec)
            phases.append(phase)
        
        # Build overlap matrix with phases
        mo = np.zeros((3, 3), dtype=complex)
        for r in range(3):
            for c in range(3):
                angle = vector_angle(axes[r], axes[c])
                overlap = np.exp(-(angle**2) / (2 * sigma**2))
                phase_factor = np.exp(1j * (phases[r] - phases[c]))
                mo[r, c] = overlap * phase_factor
        
        U, _ = qr(mo)
        mod_U = np.abs(U)
        J = float(np.imag(U[0,1] * U[1,2] * np.conj(U[0,2]) * np.conj(U[1,1])))
        
        score = float(np.linalg.norm(mod_U - V_CKM_TARGET))
        j_error = abs(J - J_TARGET)
        
        # Track best J match while keeping reasonable score
        if j_error < best_j_error and score < 0.15:
            best_j_error = j_error
            best_result = {
                'phi1': phi1, 'phi2': phi2,
                'J': J, 'score': score, 'matrix': mod_U
            }

if best_result:
    print("\nBest CP-violating configuration:")
    print("  phi1 = " + str(round(best_result['phi1'], 4)))
    print("  phi2 = " + str(round(best_result['phi2'], 4)))
    print("  Jarlskog J = " + str(best_result['J']))
    print("  Score = " + str(round(best_result['score'], 4)))
    print("  Matrix:")
    print(np.round(best_result['matrix'], 4))
else:
    print("No configuration found with J close to target")