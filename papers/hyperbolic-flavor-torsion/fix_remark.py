path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Ensure remark is declared
if r"\newtheorem{remark}" not in tex:
    tex = tex.replace(
        r"\newtheorem{conjecture}{Conjecture}",
        r"\newtheorem{conjecture}{Conjecture}" + "\n"
        r"\newtheorem{proposition}{Proposition}" + "\n"
        r"\newtheorem{remark}{Remark}"
    )
    print("Added remark theorem env.")
else:
    print("remark already declared.")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
