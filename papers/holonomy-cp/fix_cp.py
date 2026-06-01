path = r"C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp-epjc.tex"
with open(path, encoding="utf-8") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    # Remove theoremstyle lines - svjour3 handles these
    if "\\theoremstyle" in line:
        print(f"Removed: {line.rstrip()}")
        continue
    # Remove \affil lines - already replaced with \institute
    if "\\affil" in line:
        print(f"Removed: {line.rstrip()}")
        continue
    new_lines.append(line)

with open(path, "w", encoding="utf-8") as f:
    f.writelines(new_lines)
print("Done.")
