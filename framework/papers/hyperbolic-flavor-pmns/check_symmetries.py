import snappy

M = snappy.OrientableClosedCensus[0]  # m003
G = M.symmetry_group()
print("H1:", M.homology())
print("Group:", G)
print("\nEnumerating symmetries and their action on H1 (Z/5 + Z/5):\n")

for i in range(G.order()):
    g = G[i]
    print(f"Element {i}:")
    try:
        A = g.homology_matrix()
        print(A)
    except Exception as e:
        print("  No homology action:", e)
