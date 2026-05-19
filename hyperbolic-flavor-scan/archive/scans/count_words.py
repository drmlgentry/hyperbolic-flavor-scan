import snappy, itertools

M = snappy.OrientableClosedCensus[43]
G = M.fundamental_group()
base = [g for g in G.generators() if g.islower()]
all_gens = base + [g.upper() for g in base]
inv = {g: g.upper() if g.islower() else g.lower() for g in all_gens}

total = 0
for length in range(1, 11):
    count = sum(
        1 for w in itertools.product(all_gens, repeat=length)
        if all(w[i] != inv.get(w[i+1],"") for i in range(len(w)-1))
    )
    total += count
    print(f"Length {length}: {count:>8,} words  (cumulative: {total:>10,})")
