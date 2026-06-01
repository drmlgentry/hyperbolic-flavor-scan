path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()
tex = tex.replace("ΔAIC", r"$\Delta$AIC")
with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Fixed.")
