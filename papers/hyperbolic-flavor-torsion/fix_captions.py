path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Replace the two broken captions with clean versions
old1 = r"\caption{A-factor twist spectrum \$|\phi(\gamma)|\$ vs.\ word length for \mfld{m006}, coloured by homology class \\in\mathbb{Z}/5\$. The isolated class-4 point (\word{AAABAB}, \$\phi_{\rm fold}=0.0053\degr\$) is marked; the grey band is the spectral gap \$[0.005\degr,1.61\degr]\$; the dashed line marks \$\theta_{13}^{\rm CKM}=0.201\degr\$.}"

new1 = r"\caption{A-factor twist spectrum $|\phi(\gamma)|$ (folded to $[0^\circ,90^\circ]$) vs.\ word length for \mfld{m006}, coloured by homology class $k\in\mathbb{Z}/5$. The isolated class-4 point (\word{AAABAB}, $\phi_{\rm fold}=0.0053^\circ$, length~6) is annotated; the grey band marks the spectral gap $[0.005^\circ,1.61^\circ]$; the dashed line marks $\theta_{13}^{\rm CKM}=0.201^\circ$.}"

old2 = r"\caption{Spectral floors \$\phi_{\rm floor}(k)\$ per homology class \$k\in\mathbb{Z}/5\$"
new2 = r"\caption{Spectral floors $\phi_{\rm floor}(k)$ per homology class $k\in\mathbb{Z}/5$"

tex = tex.replace(old1, new1)
tex = tex.replace(old2, new2)

# Also fix any remaining \$ -> $ and \\in -> \in in captions
import re
def fix_caption(m):
    s = m.group(0)
    s = s.replace(r'\$', '$')
    s = s.replace(r'\\in', r'\in')
    s = s.replace(r'\\leq', r'\leq')
    s = s.replace(r'\\times', r'\times')
    return s

tex = re.sub(r'\\caption\{.*?\}', fix_caption, tex, flags=re.DOTALL)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Captions fixed.")
