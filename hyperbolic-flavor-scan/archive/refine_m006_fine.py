import snappy
import numpy as np
from scipy.linalg import qr, logm

# PDG CKM moduli
V_CKM_TARGET = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.0405,  0.99914]
])

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
        if norm > 1e-10:
            return vec / norm
        else:
            return np.array([1.0, 0.0, 0.0])
    except:
        return np.array([1.0, 0.0, 0.0])

def vector_angle(v1, v2):
    dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
    return np.arccos(np.abs(dot))

def get_word_homology_image(word_str, num_gens):
    """Counts net generators in the word (placeholder – will be replaced)."""
    image = np.zeros(num_gens)
    for char in word_str:
        if char.islower():
            idx = ord(char) - ord('a')
            if idx < num_gens: 
                image[idx] += 1
        elif char.isupper():
            idx = ord(char) - ord('A')
            if idx < num_gens: 
                image[idx] -= 1
    return image

def evaluate_m006(sigma, alpha=0.0, words=['aaB','aBa','AAb'], verbose=False):
    """Evaluate m006 with optional U(1) phases."""
    M = snappy.OrientableClosedCensus[43]
    rho = M.polished_holonomy()
    num_gens = len(rho.generators())
    
    axes = []
    for w in words:
        mat = rho(w)
        axes.append(matrix_to_axis_vector(mat))
    
    # Phases: if alpha>0, we need actual homology images; for now, use placeholder.
    if alpha > 0:
        # Placeholder: use random phases from earlier script (not physically motivated)
        base_phases = np.array([0.13, 0.71, 1.37, 1.93, 2.41, 3.14][:num_gens]) * alpha
        phases = []
        for w in words:
            img = get_word_homology_image(w, num_gens)
            phase = np.dot(img, base_phases)
            phases.append(phase)
    else:
        phases = [0.0, 0.0, 0.0]
    
    # Build overlap matrix
    mo = np.zeros((3,3), dtype=complex)
    for i in range(3):
        for j in range(3):
            ang = vector_angle(axes[i], axes[j])
            overlap = np.exp(-(ang**2) / (2 * sigma**2))
            phase_factor = np.exp(1j * (phases[i] - phases[j]))
            mo[i,j] = overlap * phase_factor
    
    U, _ = qr(mo)
    U_mod = np.abs(U)
    score = np.linalg.norm(U_mod - V_CKM_TARGET)
    J = float(np.imag(U[0,1] * U[1,2] * np.conj(U[0,2]) * np.conj(U[1,1])))
    
    if verbose:
        print(f"sigma={sigma:.2f}, alpha={alpha}")
        print("Mixing matrix moduli:")
        print(U_mod.round(4))
        print(f"Score = {score:.4f}, J = {J:.2e}")
    
    return score, U_mod, J

# --- Fine sigma scan ---
sigma_vals = np.arange(0.46, 0.51, 0.01)
best_score = float('inf')
best_sigma = None
best_matrix = None
best_J = None

print("Fine sigma scan for m006 with words ['aaB','aBa','AAb'] (no phases)")
print(" sigma   score      Jarlskog")
print("-" * 35)
for s in sigma_vals:
    score, U, J = evaluate_m006(s, alpha=0.0)  # no phases
    print(f"{s:.2f}    {score:.4f}    {J:.2e}")
    if score < best_score:
        best_score = score
        best_sigma = s
        best_matrix = U
        best_J = J

print("\n" + "="*50)
print(f"Best sigma: {best_sigma:.2f} with score {best_score:.4f}")
print(f"Jarlskog: {best_J:.2e}")
print("\nMixing matrix at best sigma:")
print(best_matrix.round(4))
print("\nCKM target:")
print(V_CKM_TARGET.round(4))
print("\nDifference:")
print((best_matrix - V_CKM_TARGET).round(4))

# --- If we can get abelianization, we can compute true phases ---
print("\n" + "="*50)
print("Attempting to compute abelianization for m006...")
M = snappy.OrientableClosedCensus[43]
G = M.fundamental_group()
try:
    abel = G.abelianization()
    print("abelianization result:", abel)
    # If abel is a list of integer vectors, we can use them.
    # For now, just print.
except Exception as e:
    print(f"abelianization failed: {e}")
    print("Will need to compute homology images manually from relators.")