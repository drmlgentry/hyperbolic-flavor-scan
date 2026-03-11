import snappy
import numpy as np
from itertools import product

def trace_to_length(tr):
    """Translation length from trace of SL(2,C) matrix."""
    tr = complex(tr)
    if abs(tr) <= 2:
        return 0.0   # not hyperbolic
    return 2 * np.arccosh(abs(tr)/2)

def generate_words(num_gens, max_len=4):
    """Generate all non‑empty words up to length max_len over given number of generators and their inverses.
       Returns list of tuples where each entry is an integer in [-num_gens, -1] U [1, num_gens],
       with positive numbers representing the generator (1-indexed) and negative numbers representing its inverse.
    """
    words = []
    gens = list(range(1, num_gens+1)) + list(range(-num_gens, 0))
    for length in range(1, max_len+1):
        for idx_tuple in product(gens, repeat=length):
            words.append(idx_tuple)
    return words

def word_to_matrix(word, gen_mats, inv_mats):
    """Multiply matrices in order."""
    mat = np.eye(2, dtype=complex)
    for w in word:
        if w > 0:
            mat = gen_mats[w-1] @ mat
        else:
            mat = inv_mats[-w-1] @ mat
    return mat

def trace(mat):
    return np.trace(mat)

# Locate m007
census = snappy.OrientableClosedCensus
M = None
for i in range(len(census)):
    if census[i].name() == 'm007':
        M = census[i]
        idx = i
        break

if M is None:
    print("Manifold m007 not found. Trying index 2...")
    M = census[2]   # from earlier output, index 2 is m007
    idx = 2

vol = float(M.volume())
print(f"Manifold: {M.name()} (index {idx}), volume {vol:.4f}")

# Get representation
rho = M.polished_holonomy()
gen_names = rho.generators()
num_gens = len(gen_names)
print(f"Number of generators: {num_gens}")

if num_gens < 2:
    raise ValueError("Manifold has fewer than 2 generators.")

# Use first two generators
gen_mats = [np.array(rho.SL2C(name), dtype=complex) for name in gen_names[:2]]
inv_mats = [np.linalg.inv(m) for m in gen_mats]

# Generate words up to length 4 (including all inverses)
words = generate_words(2, max_len=4)
print(f"Generated {len(words)} words up to length 4.")

traces = []
for word in words:
    try:
        mat = word_to_matrix(word, gen_mats, inv_mats)
        tr = trace(mat)
        abs_tr = abs(tr)
        if abs_tr > 2:   # keep only hyperbolic elements
            traces.append((abs_tr, word, tr))
    except Exception as e:
        continue

# Sort by |trace|
traces.sort(key=lambda x: x[0])

# Remove near duplicates
unique_traces = []
seen = set()
tol = 1e-8
for abs_tr, word, tr in traces:
    rounded = round(abs_tr, 8)
    if rounded not in seen:
        seen.add(rounded)
        unique_traces.append((abs_tr, word, tr))

print("\nThree shortest hyperbolic elements (by |trace|):")
for i in range(min(3, len(unique_traces))):
    abs_tr, word, tr = unique_traces[i]
    L = trace_to_length(tr)
    word_str = ''.join(str(w) if w>0 else f"g{-w}⁻¹" for w in word)
    print(f"Word: {word_str}, |trace|={abs_tr:.4f}, length L={L:.4f}")
    if i == 0:
        mat = word_to_matrix(word, gen_mats, inv_mats)
        print(f"Matrix:\n{mat}")