path = r"C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp-epjc.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Add LatticeFit footnote after the fitness/null test mention
# Find a good insertion point near fitness discussion
latticefit_note = (
    " The null tests were cross-validated using the open-source "
    r"\texttt{LatticeFit} package~\cite{Gentry:LatticeFit} "
    r"\cite{Gentry:LatticeFit}, which applies log-uniform null tests "
    "to arbitrary positive-valued datasets."
)

# Add bibitem
latticefit_bib = (
    "\n\\bibitem{Gentry:LatticeFit}\n"
    "M.~L. Gentry,\n"
    "``LatticeFit v0.2.0,'' Zenodo (2026).\n"
    "\\url{https://doi.org/10.5281/zenodo.19225731}\n"
)

if "LatticeFit" not in tex:
    # Add to bibliography
    tex = tex.replace(
        r"\end{thebibliography}",
        latticefit_bib + r"\end{thebibliography}"
    )
    # Add \cite near fitness value mention
    tex = tex.replace(
        "F = 0.019",
        r"F = 0.019~\cite{Gentry:LatticeFit}",
        1
    )
    print("Added LatticeFit citation")
    print("zenodo in tex:", "zenodo" in tex)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done.")
