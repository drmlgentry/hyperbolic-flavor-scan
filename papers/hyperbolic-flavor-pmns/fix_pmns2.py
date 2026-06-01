import re
path = r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Remove all \newtheorem declarations - svjour3 provides them
tex = re.sub(r'\\newtheorem\{proposition\}.*?\n', '', tex)
tex = re.sub(r'\\newtheorem\{theorem\}.*?\n', '', tex)
tex = re.sub(r'\\newtheorem\{corollary\}.*?\n', '', tex)
tex = re.sub(r'\\newtheorem\{definition\}.*?\n', '', tex)
tex = re.sub(r'\\newtheorem\{remark\}.*?\n', '', tex)
tex = re.sub(r'\\newtheorem\{lemma\}.*?\n', '', tex)

# Fix amsthm
tex = tex.replace(r"\usepackage{amsmath,amssymb,amsthm}", r"\usepackage{amsmath,amssymb}")
tex = tex.replace(r"\usepackage{amsthm}", "")

# Fix acknowledgments
tex = tex.replace(r"\begin{acknowledgments}", r"\begin{acknowledgement}")
tex = tex.replace(r"\end{acknowledgments}", r"\end{acknowledgement}")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
print("amsthm present:", "amsthm" in tex)
print("newtheorem remaining:", tex.count(r"\newtheorem"))
