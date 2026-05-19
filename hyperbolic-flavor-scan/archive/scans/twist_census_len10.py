"""
twist_census_len10.py - extend census to length 10, focus on class-3 floor stability
Key question: does the class-3 floor (1.611 deg from ABaBAb, len 6) decrease at lengths 8-10?
"""
import snappy, numpy as np, pandas as pd, itertools, os, time

DATA_DIR = r"C:\dev\hyperbolic-flavor-scan\data"
MAX_LEN  = 10
SAVE_EVERY = 500

def slc2(G, word):
    m = G.SL2C(word)
    return np.array([[complex(m[i,j]) for j in range(2)] for i in range(2)])

def get_phi_modlam(G, word):
    try:
        L  = slc2(G, word)
        ev = np.linalg.eigvals(L)
        lam = ev[np.argmax(np.abs(ev))]
        return float(np.degrees(np.angle(lam))), float(np.max(np.abs(ev)))
    except:
        return None, None

def homology_class(word, n=5):
    return sum(1 if c.islower() else -1 for c in word) % n

def gen_words(gens, max_len):
    inv = {g: g.upper() if g.islower() else g.lower() for g in gens}
    for length in range(1, max_len+1):
        for w in itertools.product(gens, repeat=length):
            if all(w[i] != inv.get(w[i+1],"") for i in range(len(w)-1)):
                yield "".join(w)

# Only run on m006 -- that's where the gap is
idx, name = 43, "m006"
M = snappy.OrientableClosedCensus[idx]
G = M.fundamental_group()
base = [g for g in G.generators() if g.islower()]
all_gens = base + [g.upper() for g in base]

print(f"{name}  vol={float(M.volume()):.4f}  H1={M.homology()}")

outfile  = os.path.join(DATA_DIR, f"twist_census_len10_{name}.csv")
ckptfile = os.path.join(DATA_DIR, f"twist_census_len10_{name}_checkpoint.txt")

# Resume from length-7 data if available
start_idx, results = 0, []
if os.path.exists(ckptfile) and os.path.exists(outfile):
    try:
        with open(ckptfile) as f: start_idx = int(f.read().strip())
        results = pd.read_csv(outfile).to_dict("records")
        print(f"Resuming from {start_idx} ({len(results)} results)")
    except: pass

words = list(gen_words(all_gens, MAX_LEN))
print(f"Total words to evaluate: {len(words):,}")
t0 = time.time()

# Track class-3 floor as we go
class3_floor = 1.611
class3_word  = "ABaBAb"

for i, word in enumerate(words):
    if i < start_idx: continue
    phi, modlam = get_phi_modlam(G, word)
    if phi is None: continue
    phi_fold = min(abs(phi)%180, 180-abs(phi)%180)
    h1 = homology_class(word)
    results.append({"word": word, "length": len(word),
                    "phi_deg": phi, "phi_fold": phi_fold,
                    "mod_lambda": modlam, "h1_class": h1})

    # Alert on class-3 floor updates
    if h1 == 3 and phi_fold < class3_floor:
        class3_floor = phi_fold
        class3_word  = word
        print(f"  *** CLASS-3 FLOOR UPDATE: {phi_fold:.6f} deg  word={word}  len={len(word)}")

    # Alert on any phi_fold < 0.5 in non-class-4
    if phi_fold < 0.5 and h1 != 4:
        print(f"  *** LOW PHI class-{h1}: {phi_fold:.4f} deg  word={word}  len={len(word)}")

    if i % SAVE_EVERY == 0 and i > start_idx:
        pd.DataFrame(results).to_csv(outfile, index=False)
        with open(ckptfile,"w") as f: f.write(str(i))
        elapsed = time.time()-t0
        rate = (i-start_idx)/max(elapsed,1)
        eta  = (len(words)-i)/rate/60
        print(f"  [{i:>8}/{len(words):>8}] {len(results)} results  "
              f"{rate:.0f} w/s  ETA {eta:.1f} min  "
              f"class-3 floor: {class3_floor:.4f} deg")

df = pd.DataFrame(results)
df.to_csv(outfile, index=False)
with open(ckptfile,"w") as f: f.write(str(len(words)))
print(f"\nComplete. {len(df)} words saved.")

print(f"\nFinal spectral floors per class:")
for k in range(5):
    sub = df[df.h1_class == k]
    floor = sub.phi_fold.min()
    word  = sub.loc[sub.phi_fold.idxmin(), "word"]
    leng  = sub.loc[sub.phi_fold.idxmin(), "length"]
    print(f"  Class [{k}]: {floor:.6f} deg  (word={word}, len={leng})")

print(f"\nGap analysis:")
floors = {k: df[df.h1_class==k].phi_fold.min() for k in range(5)}
sorted_floors = sorted(floors.values())
print(f"  Minimum: {sorted_floors[0]:.6f} deg (class {min(floors, key=floors.get)})")
print(f"  Second:  {sorted_floors[1]:.6f} deg (class {sorted(floors, key=floors.get)[1]})")
print(f"  Gap: [{sorted_floors[0]:.4f}, {sorted_floors[1]:.4f}] deg")
print(f"  theta13_CKM = 0.201 deg in gap? {sorted_floors[0] < 0.201 < sorted_floors[1]}")
print(f"  Max/min ratio: {sorted_floors[-1]/sorted_floors[0]:.0f}x")
