with open(r"C:\dev\golden_unification\papers\CORE_MASTER_v2.tex", encoding="utf-8") as f:
    tex = f.read()

# Remove all existing newtheorem lines
import re
tex = re.sub(r'\\newtheorem\{[^}]+\}[^\n]*\n', '', tex)

# Add tikz packages and clean theorem block after amsthm
old = r"\usepackage{amsmath,amssymb,amsthm}"
new = (r"\usepackage{amsmath,amssymb,amsthm}" + "\n"
       r"\usepackage{tikz}" + "\n"
       r"\usetikzlibrary{arrows.meta,positioning}" + "\n"
       r"\newtheorem{theorem}{Theorem}" + "\n"
       r"\newtheorem{proposition}[theorem]{Proposition}" + "\n"
       r"\newtheorem{lemma}[theorem]{Lemma}" + "\n"
       r"\newtheorem{corollary}[theorem]{Corollary}" + "\n"
       r"\newtheorem{definition}[theorem]{Definition}" + "\n"
       r"\newtheorem{remark}[theorem]{Remark}")
tex = tex.replace(old, new, 1)

with open(r"C:\dev\golden_unification\papers\CORE_MASTER_v2.tex", "w", encoding="utf-8") as f:
    f.write(tex)
print("Fixed. Lines:", len(tex.splitlines()))
print("theorem defined:", r"\newtheorem{theorem}" in tex)
print("tikz loaded:", r"\usepackage{tikz}" in tex)
