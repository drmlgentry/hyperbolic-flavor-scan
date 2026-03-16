import snappy
import numpy as np
from scipy.linalg import qr, logm
import warnings
warnings.filterwarnings('ignore')

M = snappy.OrientableClosedCensus[43]
print("Manifold: " + M.name())
print("Volume: " + str(round(float(M.volume()), 6)))
print("Homology: " + str(M.homology()))
print("\nKey insight: H1 = Z/5 means phases are 5th ROOTS OF UNITY!")
print("Valid phases: 0, 2π/5, 4π/5, 6π/5, 8π/5")

G = M.fundamental_group()
rho = M.polished_holonomy()

def matrix_to_axis_vector(matrix):
    try:
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
        return vec / norm if norm > 1e-10 else None
    except:
        return None

def get_word_homology_image(word_str, num_gens=2):
    image = np.zeros(num_gens)
    for char in word_str:
        if char.islower():
            idx = ord(char) - ord('a')
            if idx < num_gens: image[idx] += 1
        elif char.isupper():
            idx = ord(char) - ord('A')
            if idx < num_gens: image[idx] -= 1
    return image

def vector_angle(v1, v2):
    dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
    return np.arccos(np.abs(dot))

# Use the best words from previous scan
words = ['aaB', 'bb', 'aa']
print("\n" + "="*60)
print("Using words: " + str(words))
print("="*60)

axes = []
hom_images = []

for word in words:
    matrix = rho(word)
    tr = float(np.abs(np.trace(matrix)))
    axis = matrix_to_axis_vector(matrix)
    axes.append(axis)
    hom = get_word_homology_image(word)
    hom_images.append(hom)
    print(word + " -> hom=" + str([int(h) for h in hom]) + " |tr|=" + str(round(tr, 3)))

print("\nAxis angles:")
for i in range(3):
    for j in range(i+1, 3):
        print("  " + words[i] + " vs " + words[j] + ": " + str(round(vector_angle(axes[i], axes[j]), 4)))

# For H1 = Z/5, we have two generators [a] and [b] in the abelianization
# The U(1) character is determined by where we send each generator
# Each can be sent to any 5th root of unity: exp(2πik/5) for k=0,1,2,3,4

V_CKM_TARGET = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.0405, 0.99914]
])
J_TARGET = 3.0e-5

print("\n" + "="*60)
print("Scanning DISCRETE phases (5th roots of unity)")
print("="*60)

omega = 2 * np.pi / 5  # Fundamental phase
results = []

# k1, k2 in {0, 1, 2, 3, 4} for the two generators
for k1 in range(5):
    for k2 in range(5):
        phi1 = k1 * omega
        phi2 = k2 * omega
        char_vec = np.array([phi1, phi2])
        
        # Compute phases for each word (mod 2π effectively)
        phases = [float(np.dot(h, char_vec)) for h in hom_images]
        
        for sigma in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
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
            
            results.append({
                'k1': k1, 'k2': k2, 'sigma': sigma,
                'J': J, 'score': score, 'matrix': mod_U,
                'phases': phases
            })

print("Total discrete configurations: " + str(len(results)))

# Find non-zero J results
nonzero_J = [r for r in results if abs(r['J']) > 1e-10]
print("Configurations with |J| > 1e-10: " + str(len(nonzero_J)))

if nonzero_J:
    print("\nAll non-zero J results:")
    nonzero_J.sort(key=lambda x: -abs(x['J']))
    for r in nonzero_J[:20]:
        print("  k1=" + str(r['k1']) + " k2=" + str(r['k2']) + 
              " sig=" + str(r['sigma']) +
              " J=" + str(round(r['J'], 8)) +
              " score=" + str(round(r['score'], 4)))

# Best by score
results.sort(key=lambda x: x['score'])
best = results[0]

print("\n" + "="*60)
print("BEST by score")
print("="*60)
print("k1=" + str(best['k1']) + " k2=" + str(best['k2']) + 
      " (phases: " + str(best['k1']) + "×2π/5, " + str(best['k2']) + "×2π/5)")
print("sigma = " + str(best['sigma']))
print("J = " + str(best['J']))
print("score = " + str(round(best['score'], 4)))
print("\nMatrix:")
print(np.round(best['matrix'], 4))

# Also try adding a SMALL continuous perturbation on top of discrete
print("\n" + "="*60)
print("Adding small continuous perturbation to best discrete config")
print("="*60)

best_k1, best_k2, best_sigma = best['k1'], best['k2'], best['sigma']
base_phi1 = best_k1 * omega
base_phi2 = best_k2 * omega

perturbed_results = []
for dp1 in np.linspace(-0.3, 0.3, 30):
    for dp2 in np.linspace(-0.3, 0.3, 30):
        phi1 = base_phi1 + dp1
        phi2 = base_phi2 + dp2
        char_vec = np.array([phi1, phi2])
        phases = [float(np.dot(h, char_vec)) for h in hom_images]
        
        mo = np.zeros((3, 3), dtype=complex)
        for r in range(3):
            for c in range(3):
                angle = vector_angle(axes[r], axes[c])
                overlap = np.exp(-(angle**2) / (2 * best_sigma**2))
                phase_factor = np.exp(1j * (phases[r] - phases[c]))
                mo[r, c] = overlap * phase_factor
        
        U, _ = qr(mo)
        mod_U = np.abs(U)
        J = float(np.imag(U[0,1] * U[1,2] * np.conj(U[0,2]) * np.conj(U[1,1])))
        score = float(np.linalg.norm(mod_U - V_CKM_TARGET))
        
        perturbed_results.append({
            'dp1': dp1, 'dp2': dp2,
            'J': J, 'score': score, 'matrix': mod_U
        })

perturbed_nonzero = [r for r in perturbed_results if abs(r['J']) > 1e-6]
print("Perturbed configs with |J| > 1e-6: " + str(len(perturbed_nonzero)))

if perturbed_nonzero:
    perturbed_nonzero.sort(key=lambda x: -abs(x['J']))
    print("\nTop perturbed results by |J|:")
    for r in perturbed_nonzero[:10]:
        print("  dp1=" + str(round(r['dp1'], 3)) + 
              " dp2=" + str(round(r['dp2'], 3)) +
              " J=" + str(round(r['J'], 8)) +
              " score=" + str(round(r['score'], 4)))
    
    # Best with good score
    good_perturbed = [r for r in perturbed_nonzero if r['score'] < 0.3]
    if good_perturbed:
        good_perturbed.sort(key=lambda x: x['score'])
        best_p = good_perturbed[0]
        print("\nBest perturbed (score < 0.3):")
        print("  J = " + str(best_p['J']))
        print("  score = " + str(round(best_p['score'], 4)))
        print("  Matrix:")
        print(np.round(best_p['matrix'], 4))

perturbed_results.sort(key=lambda x: x['score'])
best_p = perturbed_results[0]
print("\nBest perturbed by score:")
print("  dp1=" + str(round(best_p['dp1'], 3)) + " dp2=" + str(round(best_p['dp2'], 3)))
print("  J = " + str(best_p['J']))
print("  score = " + str(round(best_p['score'], 4)))
print("  Matrix:")
print(np.round(best_p['matrix'], 4))
print("\nCKM target:")
print(V_CKM_TARGET)