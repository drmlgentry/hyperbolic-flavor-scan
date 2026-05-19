"""
check_lambda1.py
Cross-check: is m003 the Weeks manifold? What does the literature say about its lambda1?
Also test whether SnapPy has any direct eigenvalue computation.
"""
import snappy

# m003 identity check
M3 = snappy.OrientableClosedCensus[1]
print("m003 identification:")
print(f"  Name: {M3.name()}")
print(f"  Volume: {float(M3.volume()):.6f}")
print(f"  Homology: {M3.homology()}")
print(f"  Alexander poly: {M3.alexander_polynomial()}")

# The Weeks manifold has volume 0.942707... 
# m003 has volume 0.981369...
# These are DIFFERENT manifolds
print(f"\nWeeks manifold volume: 0.942707")
print(f"m003 volume:           {float(M3.volume()):.6f}")
print(f"Same manifold? {abs(float(M3.volume()) - 0.942707) < 0.001}")

# Check if it's the Meyerhoff manifold
print(f"\nMeyerhoff manifold volume: 0.981369")
print(f"m003 is Meyerhoff? {abs(float(M3.volume()) - 0.981369) < 0.001}")

# m006 identity
M6 = snappy.OrientableClosedCensus[43]
print(f"\nm006 identification:")
print(f"  Name: {M6.name()}")
print(f"  Volume: {float(M6.volume()):.6f}")
print(f"  Homology: {M6.homology()}")

# Try to find the actual Weeks manifold in census
print("\nSearching for Weeks manifold (vol=0.9427...) in first 50 entries:")
for i in range(50):
    try:
        M = snappy.OrientableClosedCensus[i]
        if abs(float(M.volume()) - 0.942707) < 0.001:
            print(f"  Found at index {i}: {M.name()}, vol={float(M.volume()):.6f}, H1={M.homology()}")
    except:
        break
