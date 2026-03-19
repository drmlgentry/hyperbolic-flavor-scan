"""
twist_census_len12_m006.py
Extended twist census for m006 up to word length 12.

Features:
- Reduced words only (no immediate inverse cancellation)
- Genuine loxodromics only (|λ| > 1.01)
- Homology class mod 5
- Streaming CSV output
- Running spectral floors per class
"""

import snappy
import numpy as np
import csv
import time

MAX_LEN = 12
LAM_EPS = 1.01          # filter out trivial / near-trivial
OUTFILE = r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m006.csv"

# --- Manifold ---
M = snappy.OrientableClosedCensus[43]  # m006
G = M.fundamental_group()

# Generators
letters = ["a", "b", "A", "B"]
inverse = {"a": "A", "A": "a", "b": "B", "B": "b"}


def to_numpy(m):
    return np.array([[complex(m[i, j]) for j in range(2)] for i in range(2)])


def homology_class(word, n=5):
    return sum(1 if c.islower() else -1 for c in word) % n


def phi_fold_deg(lam):
    phi = np.degrees(np.angle(lam))
    return min(abs(phi) % 180, 180 - abs(phi) % 180)


# --- Running floors per class ---
floors = {k: float("inf") for k in range(5)}

start_time = time.time()
count_total = 0
count_kept = 0

print("Starting length-12 twist census on m006...")
print("Output:", OUTFILE)

with open(OUTFILE, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["word", "length", "h1_class", "phi_fold_deg", "abs_lambda"])

    # --- Recursive generator of reduced words ---
    def extend(word, last_letter):
        global count_total, count_kept

        L = len(word)

        if L > 0:
            count_total += 1

            try:
                mat = to_numpy(G.SL2C(word))
                ev = np.linalg.eigvals(mat)
                lam = ev[np.argmax(np.abs(ev))]
                abs_lam = abs(lam)

                if abs_lam > LAM_EPS:
                    phi = phi_fold_deg(lam)
                    h1 = homology_class(word)

                    writer.writerow([word, L, h1, phi, abs_lam])
                    count_kept += 1

                    if phi < floors[h1]:
                        floors[h1] = phi
                        print(
                            f"NEW FLOOR: class {h1}  φ={phi:.6f}°  "
                            f"len={L}  word={word}"
                        )

            except Exception:
                pass

        if L == MAX_LEN:
            return

        for c in letters:
            if last_letter and c == inverse[last_letter]:
                continue
            extend(word + c, c)

    extend("", None)

elapsed = time.time() - start_time

print("\n===== DONE =====")
print(f"Total words examined: {count_total}")
print(f"Genuine loxodromics kept: {count_kept}")
print(f"Elapsed time: {elapsed:.1f} s")

print("\nFinal spectral floors (degrees):")
for k in range(5):
    val = floors[k]
    print(f"Class {k}: {val if val < float('inf') else 'none'}")