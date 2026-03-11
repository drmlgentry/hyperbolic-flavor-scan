import snappy
import numpy as np
from itertools import product

def trace_to_length(tr):
    tr = complex(tr)
    if abs(tr) <= 2:
        return 0.0
    return 2 * np.arccosh(abs(tr)/2)

def generate_words(num_gens, max_len=4):
    words = []
    gens = list(range(1, num_gens+1)) + list(range(-num_gens, 0))
    for length in range(1, max_len+1):
        for idx_tuple in product(gens, repeat=length):
            words.append(idx_tuple)
    return words

def word_to_matrix(word, gen_mats, inv_mats):
    mat = np.eye(2, dtype=complex)
    for w in word:
        if w > 0:
            mat = gen_mats[w-1] @ mat
        else:
            mat = inv_mats[-w-1] @ mat
    return mat

def trace(mat):
    return np.trace(mat)

# List of top candidate names from your scan
candidate_names = ['m007', 'm011', 'm019', 'm004']

for name in candidate_names:
    print(f"\n{'='*50}")
    print(f"Analyzing {name}")
    print('='*50)

    # Locate the manifold
    M = None
    for i in range(len(snappy.OrientableClosedCensus)):
        if snappy.OrientableClosedCensus[i].name() == name:
            M = snappy.OrientableClosedCensus[i]
            idx = i
            break

    if M is None:
        print(f"Manifold {name} not found.")
        continue

    vol = float(M.volume())
    print(f"Index: {idx}, Volume: {vol:.4f}")

    # Get representation
    rho = M.polished_holonomy()
    gen_names = rho.generators()
    num_gens = len(gen_names)
    print(f"Number of generators: {num_gens}")

    if num_gens < 2:
        print("  Skipping: fewer than 2 generators.")
        continue

    # Use first two generators (most manifolds have 2 or 3)
    gen_mats = [np.array(rho.SL2C(name), dtype=complex) for name in gen_names[:2]]
    inv_mats = [np.linalg.inv(m) for m in gen_mats]

    # Print generator traces
    print("\nGenerator traces:")
    for i, gname in enumerate(gen_names[:2]):
        tr = trace(gen_mats[i])
        L = trace_to_length(tr)
        print(f"  Generator {gname}: |trace|={abs(tr):.4f}, length L={L:.4f}")

    # Generate words up to length 4
    words = generate_words(2, max_len=4)
    print(f"Generated {len(words)} words up to length 4.")

    traces = []
    for word in words:
        try:
            mat = word_to_matrix(word, gen_mats, inv_mats)
            tr = trace(mat)
            abs_tr = abs(tr)
            if abs_tr > 2:   # keep only hyperbolic
                traces.append((abs_tr, word, tr))
        except:
            continue

    traces.sort(key=lambda x: x[0])
    # Remove near duplicates
    unique = []
    seen = set()
    for abs_tr, word, tr in traces:
        rounded = round(abs_tr, 8)
        if rounded not in seen:
            seen.add(rounded)
            unique.append((abs_tr, word, tr))

    print("\nThree shortest hyperbolic elements:")
    for i in range(min(3, len(unique))):
        abs_tr, word, tr = unique[i]
        L = trace_to_length(tr)
        word_str = ''.join(str(w) if w>0 else f"g{-w}⁻¹" for w in word)
        print(f"  Word: {word_str}, |trace|={abs_tr:.4f}, length L={L:.4f}")