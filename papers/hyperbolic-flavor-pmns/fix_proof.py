path = r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-rip.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# amsthm is already loaded but proof env needs it — check
if "amsthm" not in tex:
    tex = tex.replace(
        r"\usepackage{amsmath}",
        r"\usepackage{amsmath}" + "\n" + r"\usepackage{amsthm}"
    )
    print("Added amsthm")
else:
    print("amsthm already present — checking for conflict")
    # Ensure amsthm comes before any \newtheorem
    print("Line:", next((i+1 for i,l in enumerate(tex.splitlines()) 
                         if "amsthm" in l), "not found"))

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
