path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Replace svjour3 with article class
tex = tex.replace(
    r"\documentclass[smallextended]{svjour3}",
    r"\documentclass[12pt,a4paper]{article}"
)
tex = tex.replace(r"\smartqed", "")
tex = tex.replace(r"\journalname{Communications in Mathematical Physics}", "")

# Fix author/institute block -> standard article format
old_author = (
    r"\author{Marvin L.~Gentry}"
    "\n"
    r"\institute{Independent Researcher, Seattle, WA, USA \\"
    "\n"
    r"\email{drmlgentry@protonmail.com}}"
)
new_author = (
    r"\author{Marvin L.~Gentry\thanks{Independent Researcher, Seattle, WA, USA."
    "\n"
    r"  \texttt{drmlgentry@protonmail.com}}}"
)
tex = tex.replace(old_author, new_author, 1)

# Fix acknowledgements environment
tex = tex.replace(
    r"\begin{acknowledgements}",
    r"\section*{Acknowledgements}"
)
tex = tex.replace(r"\end{acknowledgements}", "")

# conjecture environment -> theorem-style
tex = tex.replace(r"\begin{conjecture}", 
    r"\medskip\noindent\textbf{Conjecture.}\itshape")
tex = tex.replace(r"\end{conjecture}", r"\upshape\medskip")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Switched to article class.")
