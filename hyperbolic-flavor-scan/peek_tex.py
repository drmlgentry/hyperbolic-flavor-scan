path = r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-hyperbolic-flavor-ckm.tex"
with open(path, "r", encoding="utf-8") as f:
    lines = f.readlines()

# Title, author, abstract
print("=== TITLE/AUTHOR (lines 1-50) ===")
for i, l in enumerate(lines[:50], 1):
    print(f"{i}: {l}", end="")

print("\n=== CONCLUSION (last 80 lines) ===")
for i, l in enumerate(lines[-80:], len(lines)-79):
    print(f"{i}: {l}", end="")