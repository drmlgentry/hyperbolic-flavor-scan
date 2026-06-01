import re

files = [
    r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-ckm-rip.tex",
    r"C:\dev\framework\papers\hyperbolic-flavor-pmns\gentry-pmns-rip.tex",
]

old = (r"Code: \url{https://github.com/drmlgentry/hyperbolic-flavor-scan},"
       "\n\\texttt{analysis/stronger\\_null\\_tests.py}.")

new = (r"Null tests were run using the open-source \texttt{LatticeFit} "
       r"package~\cite{Gentry:LatticeFit} (Zenodo DOI: "
       r"\href{https://doi.org/10.5281/zenodo.19225731}{10.5281/zenodo.19225731}).")

latticefit_bib = (
    "\n\\bibitem{Gentry:LatticeFit}\n"
    "M.~L.~Gentry,\n"
    "``LatticeFit v0.2.0: A Python package for detecting and validating\n"
    "discrete multiplicative structure in empirical data,''\n"
    "Zenodo (2026).\n"
    "\\href{https://doi.org/10.5281/zenodo.19225731}"
    "{doi:10.5281/zenodo.19225731}\n"
)

for path in files:
    with open(path, encoding="utf-8") as f:
        tex = f.read()
    tex = tex.replace(old, new)
    if "Gentry:LatticeFit" not in tex:
        tex = tex.replace(r"\end{thebibliography}",
                          latticefit_bib + r"\end{thebibliography}", 1)
    with open(path, "w", encoding="utf-8") as f:
        f.write(tex)
    print(f"Updated: {path.split(chr(92))[-1]}")
    print(f"  LatticeFit cite: {'Gentry:LatticeFit' in tex}")
    print(f"  Zenodo DOI: {'zenodo.19225731' in tex}")
