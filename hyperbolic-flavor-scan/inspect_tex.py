path = r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-hyperbolic-flavor-ckm.tex"
with open(path, "r", encoding="utf-8") as f:
    lines = f.readlines()

# Print lines 25-60 (abstract and kernel area)
print("=== ABSTRACT AREA (lines 25-60) ===")
for i, l in enumerate(lines[24:60], start=25):
    print(f"{i}: {repr(l)}")

print()
print("=== KERNEL AREA (lines 150-200) ===")
for i, l in enumerate(lines[149:200], start=150):
    print(f"{i}: {repr(l)}")