import snappy

M = snappy.OrientableClosedCensus[39]
G = M.symmetry_group()
isoms = G.isometries()
f = isoms[0]

# Inspect available attributes
print("Isometry attributes:")
attrs = [a for a in dir(f) if not a.startswith("__")]
for a in attrs:
    print(f"  {a}")

print()
# Try homology_matrix directly
print("Trying homology_matrix on each isometry:")
for i, g in enumerate(isoms):
    try:
        mat = g.homology_matrix()
        print(f"  Isom {i}: homology_matrix = {mat}")
        print(f"    type = {type(mat)}")
        try:
            # Try to get the integer entry
            import numpy as np
            arr = np.array(mat)
            print(f"    as array shape={arr.shape}: {arr}")
            entry = int(arr.flat[0]) % 55
            print(f"    action mod 55 = x -> {entry}x")
        except Exception as e:
            print(f"    array conversion failed: {e}")
    except Exception as e:
        print(f"  Isom {i}: homology_matrix failed: {e}")

# Also try: is_orientation_reversing or similar
print()
print("Looking for orientation-related methods:")
for a in attrs:
    if any(x in a.lower() for x in ['orient','chiral','reflect','reverse']):
        print(f"  Found: {a}")
        try:
            val = getattr(f, a)
            if callable(val):
                print(f"    callable, result: {val()}")
            else:
                print(f"    value: {val}")
        except Exception as e:
            print(f"    error: {e}")
