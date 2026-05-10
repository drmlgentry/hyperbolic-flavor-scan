path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Check which bibitems exist
import re
bibitems = re.findall(r'\\bibitem\{(\w+)\}', tex)
print("Current bibitems:", bibitems)
