import random
import numpy as np
import pandas as pd
import snappy
from pathlib import Path

# ========= CONFIG =========
N_SAMPLES = 200000
WORD_LENGTH = 80

OUT_CSV = r"C:\dev\hyperbolic-flavor-scan\data\twist_random_L80_m006.csv"

GENERATORS = ["a", "b", "A", "B"]
inv = {"a":"A","A":"a","b":"B","B":"b"}
# ==========================

M = snappy.OrientableClosedCensus[43]
G = M.fundamental_group()

rows = []

attempts = 0
loxodromic_count = 0
stored_count = 0

def random_reduced_word(L):
    w = []
    prev = None
    for _ in range(L):
        g = random.choice(GENERATORS)
        while prev and g == inv[prev]:
            g = random.choice(GENERATORS)
        w.append(g)
        prev = g
    return "".join(w)

for i in range(N_SAMPLES):

    attempts += 1
    w = random_reduced_word(WORD_LENGTH)

    try:
        Lmat = G.SL2C(w)
        mat = np.array(Lmat, dtype=np.complex128)

        if not np.all(np.isfinite(mat)):
            continue

        eig = np.linalg.eigvals(mat)
        lam = max(eig, key=lambda z: abs(z))

        if not np.isfinite(lam):
            continue

        # ---- loxodromic filter ----
        if abs(lam) <= 1.001:
            continue

        loxodromic_count += 1

        phi = np.degrees(np.angle(lam))
        phi_fold = min(abs(phi) % 180, 180 - (abs(phi) % 180))

        if not np.isfinite(phi_fold):
            continue

        # ---- homology (safe fallback) ----
        try:
            h1 = G.homology_class(w)
            h1 = int(h1)
        except:
            h1 = -1   # unknown class

        rows.append([w, WORD_LENGTH, h1, phi_fold, abs(lam)])
        stored_count += 1

    except:
        continue

    if (i + 1) % 50000 == 0:
        print(
            f"{i+1} attempts | "
            f"loxodromic: {loxodromic_count} | "
            f"stored: {stored_count}"
        )

print("\n===== DONE =====")
print(f"Attempts: {attempts}")
print(f"Loxodromic elements: {loxodromic_count}")
print(f"Stored rows: {stored_count}")

df = pd.DataFrame(
    rows,
    columns=["word", "length", "h1_class", "phi_fold_deg", "abs_lambda"]
)

Path(OUT_CSV).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_CSV, index=False)

print(f"\nSaved: {OUT_CSV}")