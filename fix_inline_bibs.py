import re

BIBITEM_CKM = r"""
\bibitem{Gentry:CKM}
M.~L. Gentry,
\newblock {CKM} quark mixing from geodesic axes of hyperbolic 3-manifold holonomy,
\newblock \emph{Phys. Rev. D} (2026), submitted March 2026, temporary {ID} es2026mar11\_966.

\bibitem{Gentry:PMNS}
M.~L. Gentry,
\newblock Lepton mixing from Borel structure of hyperbolic holonomy,
\newblock \emph{Phys. Rev. D} (2026), submitted March 2026, temporary {ID} es2026mar13\_942.
"""

inline_papers = [
    r'C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp.tex',
    r'C:\dev\framework\papers\shape-space\gentry-shape-space.tex',
]

for path in inline_papers:
    with open(path, 'r', encoding='utf-8-sig') as f:
        tex = f.read()

    if 'Gentry:CKM' in tex and 'bibitem{Gentry:CKM}' in tex:
        print(path.split('\\')[-1] + ': bibitem already present -- skipping')
        continue

    # Insert before \end{thebibliography}
    if r'\end{thebibliography}' in tex:
        tex = tex.replace(
            r'\end{thebibliography}',
            BIBITEM_CKM + r'\end{thebibliography}'
        )
        with open(path, 'w', encoding='utf-8') as f:
            f.write(tex)
        print(path.split('\\')[-1] + ': bibitem entries inserted')
    else:
        print(path.split('\\')[-1] + ': WARNING -- no thebibliography end found')

print()
print('Done. Recompile holonomy-cp and shape-space (no bibtex needed).')
print('Recompile flavor-mixing with bibtex.')
