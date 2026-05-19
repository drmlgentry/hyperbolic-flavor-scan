"""
paper9_u1_torsion.py
Test: do the 5 cosets of H1(m003)=Z/5 partition the twist spectrum into 5 bands?
For each word w, compute:
  1. phi(w) -- twist angle
  2. The homology class [w] in H1(M, Z) = Z/5
  3. Correlate: do words in the same coset have similar phi values?

The homology class of a word w = g1 g2 ... gn is:
  sum of exponents of each generator mod 5
  (since H1 = abelianization of pi1, we just count signed generator occurrences)
"""
import snappy, numpy as np, pandas as pd

def to_numpy(m):
    return np.array([[complex(m[i,j]) for j in range(2)] for i in range(2)])

def get_phi(G, word):
    try:
        L  = to_numpy(G.SL2C(word))
        ev = np.linalg.eigvals(L)
        lam = ev[np.argmax(np.abs(ev))]
        return float(np.degrees(np.angle(lam)))
    except:
        return None

def homology_class(word, n_torsion=5):
    """
    Compute homology class of word in H1 = Z/n_torsion.
    Convention: lowercase = +1, uppercase = -1 for each generator.
    Sum all exponents mod n_torsion.
    """
    total = 0
    for c in word:
        if c.islower():
            total += 1
        else:
            total -= 1
    return total % n_torsion

for idx, name in [(1,"m003"),(43,"m006")]:
    M = snappy.OrientableClosedCensus[idx]
    G = M.fundamental_group()
    print(f"\n{'='*55}")
    print(f"{name}  H1={M.homology()}")

    # Load length-7 census
    try:
        df = pd.read_csv(f"C:\\dev\\hyperbolic-flavor-scan\\data\\twist_census_len7_{name}.csv")
    except:
        print("  No length-7 data found -- run twist_census_len7.py first")
        continue

    # Deduplicate by phi_fold
    df["phi_fold"] = df["phi_deg"].apply(lambda p: min(abs(p)%180, 180-abs(p)%180))
    df["phi_key"]  = df["phi_fold"].round(3)
    df = df.sort_values(["phi_key","length","word"])
    df = df.drop_duplicates("phi_key", keep="first")

    # Compute homology class for each word
    df["h1_class"] = df["word"].apply(lambda w: homology_class(w))

    print(f"\nTwist spectrum partitioned by H1=Z/5 coset:")
    print(f"{'Coset':>8}  {'Count':>6}  {'Mean phi_fold':>14}  {'Min phi_fold':>12}  {'Max phi_fold':>12}")
    for coset in range(5):
        sub = df[df.h1_class == coset]
        if len(sub) == 0:
            print(f"  [{coset}]   empty")
            continue
        print(f"  [{coset}]   {len(sub):>4}  {sub.phi_fold.mean():>14.3f}  "
              f"{sub.phi_fold.min():>12.3f}  {sub.phi_fold.max():>12.3f}")

    print(f"\nSM targets by coset:")
    targets = [(68.0,"delta_CKM"),(33.41,"theta12_nu"),(49.1,"theta23_nu"),
               (8.54,"theta13_nu"),(2.38,"theta23_CKM"),(13.04,"theta12_CKM")]
    for phi_t, label in targets:
        near = df[abs(df.phi_fold - phi_t) < 2.0]
        if len(near):
            cosets = sorted(near.h1_class.unique())
            print(f"  {label:20s} (target={phi_t:.2f}): cosets {cosets}  "
                  f"best match phi={near.loc[near.phi_fold.sub(phi_t).abs().idxmin(),'phi_fold']:.3f}")

    print(f"\nNear-identity geodesic (phi_fold < 0.1):")
    low = df[df.phi_fold < 0.1]
    if len(low):
        print(low[["word","length","phi_fold","h1_class"]].to_string())
    else:
        print("  None")
