path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

tex = tex.replace(
    r"\section*{Acknowledgements}",
    r"\section*{Data Availability}"
    "\n"
    r"All census data and analysis scripts supporting this study are publicly"
    "\n"
    r"available at \url{https://github.com/drmlgentry/hyperbolic-flavor-scan}."
    "\n"
    r"The complete word-length-12 census files (\texttt{twist\_census\_len12\_m003.csv}"
    "\n"
    r"and \texttt{twist\_census\_len12\_m006.csv}) are included in the repository."
    "\n\n"
    r"\section*{Acknowledgements}"
)
with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Data availability statement added.")
