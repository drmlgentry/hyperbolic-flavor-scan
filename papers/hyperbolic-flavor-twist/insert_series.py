path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# 1. Add companion NPB bib entries before \end{thebibliography}
new_bibs = r"""
\bibitem{GentryNPB1}
  M.~L.~Gentry,
  ``Geometric Origin of CP Phases from Hyperbolic Holonomy,''
  Nucl.\ Phys.\ B, submitted March 2026 (NPB-S-26-00539).
\bibitem{GentryNPB2}
  M.~L.~Gentry,
  ``Discrete Mixing Operators from Boundary Sector Geometry,''
  Nucl.\ Phys.\ B, submitted March 2026 (NPB-S-26-00540).
\bibitem{GentryNPB3}
  M.~L.~Gentry,
  ``Scale-Free Quadratic Forms, Symmetric Space Geometry, and Arithmetic
  Logarithmic Lattices,''
  Nucl.\ Phys.\ B, submitted March 2026 (NPB-S-26-00538).
"""
tex = tex.replace(
    r'\end{thebibliography}',
    new_bibs + r'\end{thebibliography}'
)
print('Bib entries added')

# 2. Insert series sentence at end of first paragraph of introduction
# Find the sentence ending "...rather than generic---mixing angles."
# and append the series sentence after it.
old_intro_end = (
    r'In a companion series~\cite{GentryCKM,GentryPMNS,GentryCP} we established'
)
new_intro_start = (
    r'This paper is the fourth in a series submitted to this journal; '
    r'the companion papers~\cite{GentryNPB1,GentryNPB2,GentryNPB3} '
    r'developing the mathematical foundations are currently under review here.' + '\n\n' +
    r'In a companion series~\cite{GentryCKM,GentryPMNS,GentryCP} we established'
)
if old_intro_end in tex:
    tex = tex.replace(old_intro_end, new_intro_start, 1)
    print('Series sentence inserted in introduction')
else:
    print('WARNING: intro anchor not found -- insert manually')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done.')
