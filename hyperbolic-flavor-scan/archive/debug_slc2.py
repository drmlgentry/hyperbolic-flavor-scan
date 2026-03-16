"""
debug_slc2.py - find the correct SL2C call for closed census manifolds
"""
import snappy, numpy as np

M3 = snappy.OrientableClosedCensus[1]
print(f"Manifold: {M3.name()}")

# Method 1: fundamental_group().SL2C
G = M3.fundamental_group()
print(f"Generators string: '{G.generators()}'")
print(f"Type of G: {type(G)}")
print(f"Dir of G (SL2C-related): {[x for x in dir(G) if 'SL' in x or 'sl' in x or 'matrix' in x.lower() or 'holonomy' in x.lower()]}")

# Method 2: direct on manifold
print(f"\nDir of M (SL2C-related): {[x for x in dir(M3) if 'SL' in x or 'sl' in x or 'holonomy' in x.lower()]}")

# Method 3: the way twist_census.py does it (from the working script)
# In our working scripts we use: G = M.fundamental_group(); G.SL2C(word)
try:
    result = G.SL2C('a')
    print(f"\nG.SL2C('a') works: {result}")
except Exception as e:
    print(f"\nG.SL2C('a') failed: {e}")

# Method 4: try with word 'ab'
try:
    result = G.SL2C('ab')
    print(f"G.SL2C('ab') works: {result}")
except Exception as e:
    print(f"G.SL2C('ab') failed: {e}")

# Method 5: check if it needs snap_to_spec or similar
try:
    result = M3.holonomy_matrix('a')
    print(f"M3.holonomy_matrix('a') works: {result}")
except Exception as e:
    print(f"M3.holonomy_matrix('a') failed: {e}")

# Method 6: what twist_census.py actually does
# From our session: G = M.fundamental_group(); then G.SL2C(word)
# But OrientableClosedCensus manifolds might need: M.fundamental_group(simplify_presentation=False)
try:
    G2 = M3.fundamental_group(simplify_presentation=False)
    result = G2.SL2C('a')
    print(f"\nG(simplified=False).SL2C('a') works: {result}")
except Exception as e:
    print(f"\nG(simplified=False).SL2C('a') failed: {e}")
