path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

old = r"\acknowledgments"
new = (
    r"\subsection{Arithmetic origin}" + "\n"
    r"The Meyerhoff manifold \mfld{m003} is proven arithmetic by Chinburg et al.~\cite{Chinburg98}: "
    r"its invariant trace field has discriminant $-283$ (one complex place) and its invariant "
    r"quaternion algebra is ramified at both real places but no finite places. "
    r"The connection to $H_1=\mathbb{Z}/5$ is explicit: \mfld{m003} is derived from a "
    r"quaternion algebra \emph{as} $H_1(\mfld{m003},\mathbb{Z})=\mathbb{Z}/5$~\cite{Chinburg98}. "
    r"\mfld{m006} is expected to be arithmetic by similar criteria (small volume, torsion homology), "
    r"but its explicit invariant trace field and quaternion algebra have not been computed in the literature "
    r"to our knowledge; this is deferred to future work using the Snap package~\cite{Snap}." + "\n\n"
    r"The Ramanujan conjecture for $\mathrm{GL}(2)$ automorphic forms over the invariant trace field "
    r"would imply $\lambda_1 \geq 1$ for any congruence cover of \mfld{m003}~\cite{MaclachlanReid}. "
    r"Our estimate $\lambda_1(\mfld{m003})\approx 2.48$ is consistent with this bound." + "\n\n"
    r"\acknowledgments"
)
tex = tex.replace(old, new, 1)

# Add Chinburg reference to bibliography
old_bib = r"\bibitem{Snap}"
new_bib = (
    r"\bibitem{Chinburg98}" + "\n"
    r"  T.~Chinburg, E.~Friedman, K.~N.~Jones, A.~W.~Reid," + "\n"
    r"  ``The arithmetic hyperbolic 3-manifold of smallest volume,''" + "\n"
    r"  Ann.\ Scuola Norm.\ Sup.\ Pisa \textbf{30}, 1 (2001)." + "\n"
    r"\bibitem{Snap}"
)
tex = tex.replace(old_bib, new_bib, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Arithmetic section added.")
