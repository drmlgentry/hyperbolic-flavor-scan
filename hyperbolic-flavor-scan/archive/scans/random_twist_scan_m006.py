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
# ==========================

M = snappy.OrientableClosedCensus[43]  # m006

rows = []

def random_reduced_word(L):
    w = []
    prev = None
    inv = {"a":"A","A":"a","b":"B","B":"b"}

    for _ in range(L):
        g = random.choice(GENERATORS)
        while prev and g == inv[prev]:
            g = random.choice(GENERATORS)
        w.append(g)
        prev = g

    return "".join(w)

for i in range(N_SAMPLES):

    w = random_reduced_word(WORD_LENGTH)

    try:
        Lmat = M.fundamental_group().SL2C(w)
        eig = np.linalg.eigvals(np.array(Lmat))

        lam = max(eig, key=lambda z: abs(z))

        if abs(lam) <= 1.01:
            continue

        phi = np.degrees(np.angle(lam))
        phi_fold = min(abs(phi) % 180, 180 - (abs(phi) % 180))

        h1 = M.fundamental_group().homology_class(w)

        rows.append([w, WORD_LENGTH, h1, phi_fold, abs(lam)])

    except:
        pass

    if (i + 1) % 50000 == 0:
        print(f"{i+1} samples...")

df = pd.DataFrame(
    rows,
    columns=["word", "length", "h1_class", "phi_fold_deg", "abs_lambda"]
)

Path(OUT_CSV).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_CSV, index=False)

print(f"\nSaved: {OUT_CSV}")