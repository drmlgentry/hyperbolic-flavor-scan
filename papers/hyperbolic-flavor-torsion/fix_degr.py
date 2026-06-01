import re

path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# In figure captions, replace $...\degr...$ with $...^\circ...$
# \degr is defined as ^{\circ} which works outside math but causes issues in some contexts
# Safest fix: in any $ math $ context inside \caption{}, replace \degr with ^\circ
def fix_math_degr(m):
    return m.group(0).replace(r'\degr', r'^\circ')

# Fix all math mode occurrences of \degr -> ^\circ inside captions
tex = re.sub(r'\\caption\{.*?\}', 
             lambda m: m.group(0).replace(r'\degr', r'^\circ'),
             tex, flags=re.DOTALL)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Fixed \\degr in captions.")
