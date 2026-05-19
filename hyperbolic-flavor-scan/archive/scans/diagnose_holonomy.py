import snappy
import numpy as np

print("=== DIAGNOSING HOLONOMY EXTRACTION ===", flush=True)

M = snappy.Manifold('m003')
M.dehn_fill((-2,3))
print(f"M_PMNS vol={float(M.volume()):.5f}", flush=True)

G = M.fundamental_group(simplify_presentation=True)
print(f"Generators: {G.generators()}", flush=True)
print(f"Num relators: {len(G.relators())}", flush=True)
print(f"Relators: {G.relators()[:2]}", flush=True)

# Try different ways to get holonomy
print("\n--- Method 1: SL2C_representation ---", flush=True)
try:
    rho = G.SL2C_representation()
    print(f"rho type: {type(rho)}", flush=True)
    g = G.generators()[0]
    print(f"Generator: {g}", flush=True)
    A = rho(g)
    print(f"rho(g) type: {type(A)}", flush=True)
    print(f"rho(g) = {A}", flush=True)
    print(f"As numpy: {np.array(A, dtype=complex)}", flush=True)
except Exception as e:
    print(f"Error: {e}", flush=True)
    import traceback; traceback.print_exc()

print("\n--- Method 2: direct representation_matrix ---", flush=True)
try:
    A = M.fundamental_group().representation(1, 0)
    print(f"Type: {type(A)}", flush=True)
    print(f"Value: {A}", flush=True)
except Exception as e:
    print(f"Error: {e}", flush=True)

print("\n--- Method 3: length_spectrum + trace ---", flush=True)
try:
    spec = M.length_spectrum(3.0)
    print(f"Num geodesics: {len(spec)}", flush=True)
    item = spec[0]
    print(f"Item keys: {list(item.keys())}", flush=True)
    print(f"Item: {item}", flush=True)
    # Check if matrix is available
    for key in item.keys():
        print(f"  {key}: {item[key]}", flush=True)
except Exception as e:
    print(f"Error: {e}", flush=True)

print("\n--- Method 4: hyperbolic structure ---", flush=True)
try:
    print(f"Dir of M: {[x for x in dir(M) if 'holonomy' in x.lower() or 'repr' in x.lower() or 'matrix' in x.lower()]}", flush=True)
except Exception as e:
    print(f"Error: {e}", flush=True)

print("\n--- Method 5: fundamental_group methods ---", flush=True)
try:
    G2 = M.fundamental_group()
    print(f"Dir of G: {[x for x in dir(G2) if not x.startswith('_')]}", flush=True)
except Exception as e:
    print(f"Error: {e}", flush=True)
