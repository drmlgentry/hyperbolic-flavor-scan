import snappy

M = snappy.OrientableClosedCensus[43]
print("Name:", M.name())
print("Volume:", M.volume())
print("Homology:", M.homology())
G = M.fundamental_group()
print("Fundamental group presentation:", G)
print("Number of generators:", G.num_generators())
print("Relators:", G.relators())
print("Is arithmetic?", M.is_arithmetic())