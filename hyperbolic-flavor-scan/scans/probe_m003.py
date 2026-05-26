import snappy

print("SnapPy version:", snappy.__version__)

M = snappy.Manifold("m003(-2,3)")

print("\n=== Basic invariants ===")
print("Volume:", M.volume())
print("Homology:", M.homology())
print("Fundamental group:")
print(M.fundamental_group())

print("\n=== Available arithmetic methods ===")
for name in dir(M):
    if "trace" in name.lower() or "field" in name.lower() or "arith" in name.lower():
        print(name)
