path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

old_bib = r"\bibitem{MaclachlanReid}"
new_bib = (
    r"\bibitem{BalogBiggs23}" + "\n"
    r"  A.~Balog, A.~Biggs, A.~de~Faveri, A.~Harper, K.~Pushkarev, S.~Steiner," + "\n"
    r"  ``Ambient Prime Geodesic Theorems on Hyperbolic 3-Manifolds,''" + "\n"
    r"  Int.\ Math.\ Res.\ Not.\ \textbf{2023}, no.~1, 588--606." + "\n"
    r"\bibitem{Margulis69}" + "\n"
    r"  G.~A.~Margulis," + "\n"
    r"  ``Certain applications of ergodic theory to the investigation of manifolds" + "\n"
    r"  of negative curvature,'' Funct.\ Anal.\ Appl.\ \textbf{3}, 335--336 (1969)." + "\n"
    r"\bibitem{Sunada85}" + "\n"
    r"  T.~Sunada," + "\n"
    r"  ``Riemannian coverings and isospectral manifolds,''" + "\n"
    r"  Ann.\ Math.\ \textbf{121}, 169--186 (1985)." + "\n"
    r"\bibitem{MaclachlanReid}"
)
tex = tex.replace(old_bib, new_bib, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("References added.")
