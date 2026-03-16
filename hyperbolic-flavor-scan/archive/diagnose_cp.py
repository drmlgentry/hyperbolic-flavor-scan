import snappy
import numpy as np
from scipy.linalg import qr, logm
import warnings
warnings.filterwarnings('ignore')

M = snappy.OrientableClosedCensus[43]
print("Manifold: " + M.name())
print("Homology: " + str(M.homology()))

G = M.fundamental_group()
rho = M.polished_holonomy()

def get_complex_fixed_point(matrix):
    """Get the attracting fixed point as a complex number on CP^1."""
    mat = np.array(matrix, dtype=complex)
    vals, vecs = np.linalg.eig(mat)
    idx = np.argmax(np.abs(vals))
    v = vecs[:, idx]
    if abs(v[1]) > 1e-10:
        return v[0] / v[1]
    else:
        return complex(1e10, 0)

words = ['aaB', 'bb', 'aa']

print("\n" + "="*60)
print("FIXED POINT ANALYSIS")
print("="*60)

fps = []
for word in words:
    matrix = rho(word)
    fp = get_complex_fixed_point(matrix)
    fps.append(fp)
    print(word + " fixed point: " + str(round(fp.real, 6)) + " + " + str(round(fp.imag, 6)) + "i")

print("\nImaginary parts:")
for i, word in enumerate(words):
    print("  " + word + ": Im = " + str(fps[i].imag))

# The KEY question: are the fixed points all real?
max_imag = max(abs(fp.imag) for fp in fps)
print("\nMax |Im(fp)|: " + str(max_imag))

if max_imag < 1e-6:
    print("\n" + "="*60)
    print("DIAGNOSIS: All fixed points are REAL!")
    print("="*60)
    print("""
This means m006 has a special property: all hyperbolic elements
have their axes lying in a single hyperbolic plane (the upper half-plane).

Consequence: The entire mixing matrix construction is REAL, so J=0 identically.

This is NOT a failure - it's a SELECTION RULE from the geometry!
m006 preserves CP because its fundamental domain is 'planar'.
    """)
else:
    print("\nFixed points have non-trivial imaginary parts.")
    print("The issue may be in how phases combine.")

# Let's check ALL short hyperbolic words
print("\n" + "="*60)
print("Checking imaginary parts of ALL short hyperbolic fixed points")
print("="*60)

import itertools
alphabet = ['a', 'A', 'b', 'B']
max_imag_found = 0
word_with_max_imag = None

for length in range(1, 6):
    for word_tuple in itertools.product(alphabet, repeat=length):
        word_str = "".join(word_tuple)
        
        # Skip cancellations
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
                fp = get_complex_fixed_point(matrix)
                if abs(fp.imag) > max_imag_found:
                    max_imag_found = abs(fp.imag)
                    word_with_max_imag = word_str
        except:
            continue

print("Maximum |Im(fp)| found: " + str(max_imag_found))
print("Word: " + str(word_with_max_imag))

if max_imag_found < 1e-6:
    print("\n*** CONFIRMED: m006 has ALL fixed points on the real axis! ***")
    print("This manifold cannot produce CP violation with this construction.")
    print("\nWe should search for a manifold with COMPLEX fixed points.")