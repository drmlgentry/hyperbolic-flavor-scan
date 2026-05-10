path = r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-ckm-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()
tex = tex.replace(
    r"\usepackage{amsmath,amssymb,amsthm}",
    r"\usepackage{amsmath,amssymb}"
)
with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.", "amsthm" in tex)
