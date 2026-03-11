import snappy
import numpy as np
import pandas as pd
from scipy.linalg import qr
import itertools
import re
import sys

# --- CONFIGURATION & TARGETS ---
V_CKM_TARGET = np.array([
    [0.97427, 0.22536, 0.00355],
    [0.22522, 0.97339, 0.04108],
    [0.00886, 0.0405, 0.99914]
])
J_TARGET = 3.0e-5

def attracting_fixed_point(matrix):
    try:
        eigenvals, eigenvecs = np.linalg.eig(matrix)
        idx = np.argmax(np.abs(eigenvals))
        v = eigenvecs[:, idx]
        if np.abs(v[1]) < 1e-12:
            return complex(1e10, 0)
        return v[0] / v[1]
    except:
        return complex(0, 0)

def parse_homology(M):
    """Extracts free rank and torsion from M.homology() string."""
    hom_str = str(M.homology())
    free_rank = hom_str.count('Z') - hom_str.count('Z/')
    torsion = [int(t) for t in re.findall(r'Z/(\d+)', hom_str)]
    return free_rank, torsion

def get_word_homology_image(word_str, num_gens):
    """Maps a word to generator counts (abelianization image)."""
    image = np.zeros(num_gens)
    for char in word_str:
        if char.islower():
            idx = ord(char) - ord('a')
            if idx < num_gens: image[idx] += 1
        elif char.isupper():
            idx = ord(char) - ord('A')
            if idx < num_gens: image[idx] -= 1
    return image

def get_shortest_hyperbolic_elements(M, max_len=5, num_elements=3):
    try:
        G = M.fundamental_group()
        rho = M.polished_holonomy()
        num_gens = len(G.generators())
    except:
        return []
    
    alphabet = [chr(97 + i) for i in range(num_gens)] + [chr(65 + i) for i in range(num_gens)]
    hyperbolic_elements = []
    
    for length in range(1, max_len + 1):
        for word_tuple in itertools.product(alphabet, repeat=length):
            word_str = "".join(word_tuple)
            # Skip immediate cancellations
            if any(word_str[i].lower() == word_str[i+1].lower() and word_str[i] != word_str[i+1] 
                   for i in range(len(word_str)-1)):
                continue
            
            try:
                matrix = rho(word_str)
                tr = np.abs(np.trace(matrix))
                if tr > 2.001:
                    hyperbolic_elements.append({'word': word_str, 'trace': tr, 'matrix': matrix})
            except:
                continue
        if len(hyperbolic_elements) >= num_elements * 5: break

    hyperbolic_elements.sort(key=lambda x: x['trace'])
    return hyperbolic_elements[:num_elements]

def run_geometric_scan(num_manifolds=50, sigma=0.3):
    results = []
    print(f"Starting Scan: sigma={sigma}, N={num_manifolds}")
    print("-" * 60)

    for i in range(num_manifolds):
        try:
            M = snappy.OrientableClosedCensus[i]
            elements = get_shortest_hyperbolic_elements(M, max_len=6, num_elements=3)
            if len(elements) < 3: continue

            num_gens = len(M.fundamental_group().generators())
            _, torsion_list = parse_homology(M)
            
            # Phase vector for CP violation
            char_vec = np.zeros(num_gens)
            for j in range(num_gens):
                if j < len(torsion_list):
                    char_vec[j] = 2 * np.pi * np.random.randint(0, torsion_list[j]) / torsion_list[j]
                else:
                    char_vec[j] = np.random.uniform(0, 2 * np.pi)

            # Build Overlap Matrix
            M_overlap = np.zeros((3, 3), dtype=complex)
            fps = [attracting_fixed_point(el['matrix']) for el in elements]
            phases = [np.dot(get_word_homology_image(el['word'], num_gens), char_vec) for el in elements]

            for r in range(3):
                for c in range(3):
                    z1, z2 = fps[r], fps[c]
                    dist = 2.0 * np.abs(z1 - z2) / (np.sqrt(1 + np.abs(z1)**2) * np.sqrt(1 + np.abs(z2)**2))
                    M_overlap[r, c] = np.exp(-(dist**2)/(2*sigma**2)) * np.exp(1j*(phases[r]-phases[c]))

            U, _ = qr(M_overlap)
            mod_U = np.abs(U)
            J = np.imag(U[0,1]*U[1,2]*np.conj(U[0,2])*np.conj(U[1,1]))
            
            score = np.linalg.norm(mod_U - V_CKM_TARGET) + 0.1 * (np.abs(J - J_TARGET) / J_TARGET)

            results.append({
                'name': M.name(), 'vol': M.volume(), 'score': score, 'j_inv': J,
                'u11': mod_U[0,0], 'u12': mod_U[0,1], 'u13': mod_U[0,2],
                'u21': mod_U[1,0], 'u22': mod_U[1,1], 'u23': mod_U[1,2],
                'u31': mod_U[2,0], 'u32': mod_U[2,1], 'u33': mod_U[2,2]
            })
            print(f"{M.name()}: Score={score:.4f}, J={J:.2e}")
        except:
            continue

    if not results:
        print("No valid results found.")
        return
    
    df = pd.DataFrame(results).sort_values('score')
    df.to_csv("scan_results_cp_final.csv", index=False)
    print("\nTOP CANDIDATE:")
    print(df.iloc[0][['name', 'score', 'j_inv']])

def analyze_manifold(name, sigma=0.3, trials=50):
    M = snappy.OrientableClosedCensus[name]
    elements = get_shortest_hyperbolic_elements(M, max_len=6, num_elements=3)
    num_gens = len(M.fundamental_group().generators())
    _, torsion = parse_homology(M)
    
    print(f"Manifold: {name}, Homology: {M.homology()}")
    traces = [f"{el['trace']:.4f}" for el in elements]
    print(f"Words: {[el['word'] for el in elements]}, Traces: {traces}")

    best_s = 1e10
    for _ in range(trials):
        cv = np.array([np.random.uniform(0, 2*np.pi) for _ in range(num_gens)])
        fps = [attracting_fixed_point(el['matrix']) for el in elements]
        ps = [np.dot(get_word_homology_image(el['word'], num_gens), cv) for el in elements]
        mo = np.zeros((3,3), dtype=complex)
        for r in range(3):
            for c in range(3):
                d = 2.0*np.abs(fps[r]-fps[c])/(np.sqrt(1+np.abs(fps[r])**2)*np.sqrt(1+np.abs(fps[c])**2))
                mo[r,c] = np.exp(-d**2/(2*sigma**2)) * np.exp(1j*(ps[r]-ps[c]))
        U, _ = qr(mo)
        J = np.imag(U[0,1]*U[1,2]*np.conj(U[0,2])*np.conj(U[1,1]))
        s = np.linalg.norm(np.abs(U)-V_CKM_TARGET) + 0.1 * np.abs(J-J_TARGET)/J_TARGET
        if s < best_s: best_s, best_U, best_J = s, np.abs(U), J
    
    print(f"Best Score: {best_s:.4f}, J: {best_J:.4e}")
    print("Moduli Matrix:\n", best_U)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "analyze":
        name = sys.argv[2] if len(sys.argv) > 2 else "m007"
        analyze_manifold(name)
    else:
        run_geometric_scan()