import re
path = r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-ckm-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Fix 1: remove amsthm (svjour3 provides its own)
tex = tex.replace(r"\usepackage{amsthm}", "")
tex = tex.replace("% amsthm removed - svjour3 provides proof env", "")

# Fix 2: acknowledgments -> acknowledgement (svjour3 spelling)
tex = tex.replace(r"\begin{acknowledgments}", r"\begin{acknowledgement}")
tex = tex.replace(r"\end{acknowledgments}", r"\end{acknowledgement}")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
print("amsthm present:", r"\usepackage{amsthm}" in tex)
print("acknowledgment fixed:", "acknowledgement" in tex)
