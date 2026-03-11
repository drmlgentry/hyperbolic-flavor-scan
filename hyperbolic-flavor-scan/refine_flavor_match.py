import snappy
import numpy as np
from scipy.linalg import qr, logm
import itertools

# --- TARGETS ---
V_CKM_TARGET = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.0405, 0.99914]
])
J_TARGET = 3.0e-5

def get_word_homology_image(word_str, num_gens):
    """Abelianization: counts net generators in the word."""
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

def matrix_to_axis_vector(matrix):
    """
    Extract a direction vector from an SL(2,C) matrix.
    Uses the matrix logarithm to get the "axis" direction.
    """
    try:
        # Convert to numpy array with complex entries
        mat = np.array(matrix, dtype=complex)
        
        # Normalize to SL(2,C)
        det = np.linalg.det(mat)
        mat = mat / np.sqrt(det)
        
        # Take matrix logarithm
        log_mat = logm(mat)
        
        # The "direction" is encoded in the off-diagonal structure
        # Extract a characteristic vector
        a = log_mat[0, 0]
        b = log_mat[0, 1]
        c = log_mat[1, 0]
        d = log_mat[1, 1]
        
        # Create a 3-vector from the traceless part
        # Using the Pauli decomposition: M = t*I + x*σ1 + y*σ2 + z*σ3
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
    """Angle between two unit vectors in radians."""
    dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
    return np.arccos(np.abs(dot))  # Use abs to treat parallel and antiparallel as same

def find_distinct_hyperbolic_words(M, max_length=7, num_elements=3, min_angle=0.3, debug=False):
    """
    Finds hyperbolic words with geometrically distinct axes.
    """
    try:
        G = M.fundamental_group()
        rho = M.polished_holonomy()
        gens = G.generators()
        num_gens = len(gens)
    except Exception as e:
        if debug: 
            print("   Error: " + str(e))
        return [], 0
    
    if debug:
        print("   Generators: " + str(gens))
    
    alphabet = []
    for i in range(num_gens):
        alphabet.append(chr(97 + i))
        alphabet.append(chr(65 + i))
    
    candidates = []
    
    for length in range(1, max_length + 1):
        for word_tuple in itertools.product(alphabet, repeat=length):
            word_str = "".join(word_tuple)
            
            # Skip immediate cancellations
            has_cancel = False
            for i in range(len(word_str) - 1):
                c1, c2 = word_str[i], word_str[i+1]
                if c1.lower() == c2.lower() and c1 != c2:
                    has_cancel = True
                    break
            if has_cancel:
                continue
            
            try:
                matrix = rho(word_str)
                tr = float(np.abs(np.trace(matrix)))
                
                if tr > 2.01:
                    axis = matrix_to_axis_vector(matrix)
                    candidates.append({
                        'word': word_str,
                        'matrix': matrix,
                        'trace': tr,
                        'axis': axis
                    })
            except:
                continue
        
        if len(candidates) >= 100:
            break
    
    if debug:
        print("   Found " + str(len(candidates)) + " hyperbolic candidates")
    
    # Sort by trace (length proxy)
    candidates.sort(key=lambda x: x['trace'])
    
    # Greedily select elements with distinct axes
    selected = []
    for cand in candidates:
        is_distinct = True
        for sel in selected:
            angle = vector_angle(cand['axis'], sel['axis'])
            if angle < min_angle:
                is_distinct = False
                break
        
        if is_distinct:
            selected.append(cand)
            if debug:
                ax = cand['axis']
                print("   Selected: '" + cand['word'] + "' |tr|=" + str(round(cand['trace'], 4)) + 
                      " axis=[" + str(round(ax[0], 3)) + "," + 
                      str(round(ax[1], 3)) + "," + str(round(ax[2], 3)) + "]")
        
        if len(selected) >= num_elements:
            break
    
    return selected, num_gens

def evaluate_manifold(M, sigma, alpha, min_angle=0.3, debug=False):
    """
    Evaluates a manifold using axis-based geometry.
    """
    name = M.name()
    
    elements, num_gens = find_distinct_hyperbolic_words(
        M, max_length=7, num_elements=3, min_angle=min_angle, debug=debug
    )
    
    if len(elements) < 3:
        if debug:
            print("   Only found " + str(len(elements)) + " distinct elements")
        return None
    
    el = elements[:3]
    
    # Use axis vectors to compute overlaps
    axes = [e['axis'] for e in el]
    
    if debug:
        print("   Axis angles (radians):")
        for i in range(3):
            for j in range(i+1, 3):
                ang = vector_angle(axes[i], axes[j])
                print("     angle(" + str(i) + "," + str(j) + ") = " + str(round(ang, 4)))
    
    # Compute phases from homology
    char_vec = np.array([0.13, 0.71, 1.37, 1.93, 2.41, 3.14][:num_gens]) * alpha
    phases = []
    for e in el:
        hom_img = get_word_homology_image(e['word'], num_gens)
        phase = np.dot(hom_img, char_vec)
        phases.append(phase)
    
    # Build overlap matrix using axis dot products
    mo = np.zeros((3, 3), dtype=complex)
    for r in range(3):
        for c in range(3):
            # Use angle between axes
            angle = vector_angle(axes[r], axes[c])
            # Gaussian overlap based on angle
            overlap = np.exp(-(angle**2) / (2 * sigma**2))
            phase_factor = np.exp(1j * (phases[r] - phases[c]))
            mo[r, c] = overlap * phase_factor
    
    if debug:
        print("   Overlap matrix moduli:")
        print(np.round(np.abs(mo), 4))
    
    # QR orthogonalization
    U, _ = qr(mo)
    mod_U = np.abs(U)
    
    # Jarlskog invariant
    J = float(np.imag(U[0,1] * U[1,2] * np.conj(U[0,2]) * np.conj(U[1,1])))
    
    # Score = Frobenius distance to CKM
    score = float(np.linalg.norm(mod_U - V_CKM_TARGET))
    
    return score, J, mod_U, [e['word'] for e in el]

# --- MAIN ---
print("=" * 70)
print("HYPERBOLIC FLAVOR SCAN - Axis-Based Geometry")
print("=" * 70)

# First, let's examine a few manifolds in detail
print("\nStep 1: Detailed manifold analysis")
print("-" * 70)

for idx in [0, 2, 5, 10]:
    M = snappy.OrientableClosedCensus[idx]
    print("\nManifold index " + str(idx) + ": " + M.name())
    print("Volume: " + str(round(float(M.volume()), 6)))
    print("Homology: " + str(M.homology()))
    
    result = evaluate_manifold(M, sigma=0.5, alpha=0.1, min_angle=0.2, debug=True)
    if result:
        score, J, mod_U, words = result
        print("Result: Score=" + str(round(score, 4)) + ", J=" + str(round(J, 6)))
        print("Mixing matrix:")
        print(np.round(mod_U, 4))
    else:
        print("Could not find 3 distinct geodesics")

# Parameter sweep
print("\n" + "=" * 70)
print("Step 2: Parameter Sweep (first 50 manifolds)")
print("=" * 70)

results = []

for idx in range(50):
    M = snappy.OrientableClosedCensus[idx]
    name = M.name()
    
    for sig in [0.3, 0.5, 0.8, 1.2]:
        for alp in [0.1, 0.05, 0.01]:
            for ang in [0.15, 0.25, 0.4]:
                result = evaluate_manifold(M, sigma=sig, alpha=alp, min_angle=ang, debug=False)
                if result:
                    score, J, mod_U, words = result
                    results.append({
                        'idx': idx,
                        'name': name,
                        'sigma': sig,
                        'alpha': alp,
                        'min_angle': ang,
                        'score': score,
                        'J': J,
                        'matrix': mod_U,
                        'words': words
                    })

print("Total valid configurations: " + str(len(results)))

if len(results) > 0:
    # Sort by score
    results.sort(key=lambda x: x['score'])
    
    # Print top 20
    print("\nTop 20 results:")
    print("-" * 80)
    header = "Rank | Idx | Name     | sigma | alpha | ang  | Score  | Jarlskog"
    print(header)
    print("-" * 80)
    
    for i, r in enumerate(results[:20]):
        row = (str(i+1).ljust(4) + " | " +
               str(r['idx']).ljust(3) + " | " +
               r['name'].ljust(8) + " | " +
               str(r['sigma']).ljust(5) + " | " +
               str(r['alpha']).ljust(5) + " | " +
               str(r['min_angle']).ljust(4) + " | " +
               str(round(r['score'], 4)).ljust(6) + " | " +
               str(round(r['J'], 6)))
        print(row)
    
    # Best result details
    best = results[0]
    print("\n" + "=" * 70)
    print("BEST RESULT")
    print("=" * 70)
    print("Manifold: " + best['name'] + " (index " + str(best['idx']) + ")")
    print("Parameters: sigma=" + str(best['sigma']) + 
          ", alpha=" + str(best['alpha']) + 
          ", min_angle=" + str(best['min_angle']))
    print("Words: " + str(best['words']))
    print("Score: " + str(round(best['score'], 4)))
    print("Jarlskog: " + str(best['J']) + " (target: " + str(J_TARGET) + ")")
    print("\nMixing matrix |U|:")
    print(np.round(best['matrix'], 4))
    print("\nCKM target:")
    print(V_CKM_TARGET)
    print("\nDifference (|U| - CKM):")
    print(np.round(best['matrix'] - V_CKM_TARGET, 4))
    
    # Check if any result has J close to target
    print("\n" + "-" * 70)
    print("Results with Jarlskog closest to target (" + str(J_TARGET) + "):")
    j_sorted = sorted(results, key=lambda x: abs(x['J'] - J_TARGET))
    for r in j_sorted[:5]:
        print("  " + r['name'] + " idx=" + str(r['idx']) + 
              " J=" + str(round(r['J'], 8)) + 
              " score=" + str(round(r['score'], 4)))
else:
    print("\nNo valid results found!")