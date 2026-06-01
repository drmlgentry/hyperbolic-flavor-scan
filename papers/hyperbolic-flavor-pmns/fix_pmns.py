import re
path = r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Fix amsthm conflict
tex = tex.replace(r"\usepackage{amsmath,amssymb,amsthm}", r"\usepackage{amsmath,amssymb}")
tex = tex.replace(r"\usepackage{amsthm}", "")

# Fix acknowledgments spelling
tex = tex.replace(r"\begin{acknowledgments}", r"\begin{acknowledgement}")
tex = tex.replace(r"\end{acknowledgments}", r"\end{acknowledgement}")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
print("amsthm present:", "amsthm" in tex)
