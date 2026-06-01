path = r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-ckm-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()
# svjour3 already defines \proof - remove amsthm redefinition
tex = tex.replace(r"\usepackage{amsthm}", "% amsthm removed - svjour3 provides proof env")
# Also remove any \newtheorem for proof if present
import re
tex = re.sub(r'\\newtheorem\{proof\}.*?\n', '', tex)
with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
