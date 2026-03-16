"""
twist_census_len7_v4.py - correct SimpleMatrix conversion
"""
import snappy, numpy as np, pandas as pd, itertools, os, time

MANIFOLDS = [(1,"m003"),(43,"m006")]
MAX_LEN   = 7
SAVE_EVERY = 200
DATA_DIR   = r"C:\dev\hyperbolic-flavor-scan\data"

def slc2_to_numpy(mat):
    return np.array([[complex(mat[i,j]) for j in range(2)] for i in range(2)])

def get_phi_modlam(G, word):
    try:
        L  = slc2_to_numpy(G.SL2C(word))
        ev = np.linalg.eigvals(L)
        lam = ev[np.argmax(np.abs(ev))]
        return float(np.degrees(np.angle(lam))), float(np.max(np.abs(ev)))
    except:
        return None, None

def gen_words(gens, max_len):
    inv = {g: g.upper() if g.islower() else g.lower() for g in gens}
    for length in range(1, max_len+1):
        for w in itertools.product(gens, repeat=length):
            ok = True
            for i in range(len(w)-1):
                if w[i] == inv.get(w[i+1],"___"):
                    ok = False; break
            if ok:
                yield "".join(w)

for idx, name in MANIFOLDS:
    M = snappy.OrientableClosedCensus[idx]
    G = M.fundamental_group()
    base = [g for g in G.generators() if g.islower()]
    all_gens = base + [g.upper() for g in base]

    print(f"\n{'='*55}")
    print(f"{name}  vol={float(M.volume()):.4f}  H1={M.homology()}")
    print(f"Generators: {all_gens}")

    phi_test, _ = get_phi_modlam(G, all_gens[0])
    if phi_test is None:
        print("ERROR: SL2C failing -- skipping"); continue
    print(f"SL2C check: phi({all_gens[0]}) = {phi_test:.3f} deg -- OK")

    outfile   = os.path.join(DATA_DIR, f"twist_census_len7_{name}.csv")
    ckpt_file = os.path.join(DATA_DIR, f"twist_census_len7_{name}_checkpoint.txt")

    start_idx, results = 0, []
    if os.path.exists(ckpt_file) and os.path.exists(outfile):
        try:
            with open(ckpt_file) as f: start_idx = int(f.read().strip())
            results = pd.read_csv(outfile).to_dict("records")
            print(f"Resuming from {start_idx} ({len(results)} results)")
        except: start_idx, results = 0, []

    words = list(gen_words(all_gens, MAX_LEN))
    print(f"Total words: {len(words):,}")
    t0 = time.time()

    for i, word in enumerate(words):
        if i < start_idx: continue
        phi, modlam = get_phi_modlam(G, word)
        if phi is None: continue
        phi_fold = min(abs(phi)%180, 180-abs(phi)%180)
        results.append({"word": word, "length": len(word),
                        "phi_deg": phi, "phi_fold": phi_fold,
                        "mod_lambda": modlam})

        if phi_fold < 0.5:
            print(f"  *** phi_fold={phi_fold:.4f} deg  word={word}  len={len(word)}")

        if i % SAVE_EVERY == 0 and i > start_idx:
            pd.DataFrame(results).to_csv(outfile, index=False)
            with open(ckpt_file,"w") as f: f.write(str(i))
            elapsed = time.time()-t0
            rate = (i-start_idx)/max(elapsed,1)
            eta  = (len(words)-i)/rate/3600
            print(f"  [{i:>7}/{len(words):>7}] {len(results)} results "
                  f"{rate:.0f} w/s  ETA {eta:.1f}h")

    if not results:
        print("No results."); continue

    df = pd.DataFrame(results)
    df.to_csv(outfile, index=False)
    with open(ckpt_file,"w") as f: f.write(str(len(words)))
    print(f"\nComplete. {len(df)} words saved to {outfile}")

    print(f"\nSmallest phi_fold (top 10):")
    print(df.nsmallest(10,"phi_fold")[
        ["word","length","phi_deg","phi_fold"]].to_string())

    low = df[df.phi_fold < 1.0].sort_values("phi_fold")
    print(f"\nAll phi_fold < 1.0 deg ({len(low)} found):")
    if len(low):
        print(low[["word","length","phi_deg","phi_fold","mod_lambda"]].to_string())
    else:
        print("  None at length <= 7.")

print("\nDone.")
