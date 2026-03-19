"""
twist_census_len12_m003.py - length-12 twist census for Meyerhoff manifold m003
"""
import snappy, numpy as np, pandas as pd, itertools, os, time

M = snappy.OrientableClosedCensus[1]   # m003 = Meyerhoff
G = M.fundamental_group()
base = [g for g in G.generators() if g.islower()]
all_gens = base + [g.upper() for g in base]
inv = {g: g.upper() if g.islower() else g.lower() for g in all_gens}

print(f"m003 Meyerhoff  vol={float(M.volume()):.4f}  H1={M.homology()}")

def slc2(word):
    m = G.SL2C(word)
    return np.array([[complex(m[i,j]) for j in range(2)] for i in range(2)])

def homology_class(word, n=5):
    return sum(1 if c.islower() else -1 for c in word) % n

outfile = r"C:\dev\hyperbolic-flavor-scan\data\twist_census_len12_m003.csv"
results = []
t0 = time.time()
count = 0
MIN_FLOOR = {k: 999.0 for k in range(5)}

def gen_words(max_len):
    for length in range(1, max_len+1):
        for w in itertools.product(all_gens, repeat=length):
            if all(w[i] != inv.get(w[i+1],"") for i in range(len(w)-1)):
                yield "".join(w)

for word in gen_words(12):
    count += 1
    try:
        L = slc2(word)
        ev = np.linalg.eigvals(L)
        lam = ev[np.argmax(np.abs(ev))]
        al = float(np.max(np.abs(ev)))
        if al <= 1.01: continue  # skip trivial
        phi = float(np.degrees(np.angle(lam)))
        phi_fold = min(abs(phi)%180, 180-abs(phi)%180)
        h1 = homology_class(word)
        results.append({"word": word, "length": len(word), "h1_class": h1,
                        "phi_fold_deg": phi_fold, "abs_lambda": al})
        if phi_fold < MIN_FLOOR[h1]:
            MIN_FLOOR[h1] = phi_fold
            print(f"NEW FLOOR: class {h1}  phi={phi_fold:.6f}  len={len(word)}  word={word}")
    except: pass
    if count % 50000 == 0:
        elapsed = time.time()-t0
        print(f"[{count:>8}]  {elapsed:.1f}s  floors: {[f'{MIN_FLOOR[k]:.4f}' for k in range(5)]}")

df = pd.DataFrame(results)
df.to_csv(outfile, index=False)
print(f"\n===== DONE =====")
print(f"Total words: {count}  Genuine loxodromics: {len(df)}  Time: {time.time()-t0:.1f}s")
print("\nFinal spectral floors:")
for k in range(5):
    print(f"Class {k}: {MIN_FLOOR[k]:.8f}")
