import snappy, numpy as np, ast

M = snappy.OrientableClosedCensus[1]
G = M.fundamental_group()

# Show exactly what generators() returns
raw = G.generators()
print(f"raw type: {type(raw)}")
print(f"raw repr: {repr(raw)}")

# Parse it
try:
    parsed = ast.literal_eval(raw)
    print(f"ast parsed: {parsed}")
    base = [g for g in parsed if len(g)==1 and g.islower()]
except Exception as e:
    print(f"ast failed: {e}")
    base = [g for g in str(raw) if g.islower() and g.isalpha()]

all_gens = base + [g.upper() for g in base]
print(f"all_gens: {all_gens}")
print(f"first gen: {repr(all_gens[0])}")

# Try SL2C directly -- no try/except
L = np.array(G.SL2C(all_gens[0]), dtype=complex)
ev = np.linalg.eigvals(L)
lam = ev[np.argmax(np.abs(ev))]
phi = float(np.degrees(np.angle(lam)))
print(f"phi({all_gens[0]}) = {phi:.4f} deg -- SUCCESS")
