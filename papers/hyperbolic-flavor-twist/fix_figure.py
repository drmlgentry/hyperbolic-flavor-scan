path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Insert figure environment after \label{fig:spectrum} reference in text
# Find the caption reference and insert the full figure block before the next subsection
old_anchor = r'\subsection{The \texorpdfstring{$\bar{m}_b/\bar{m}_c$}{mb/mc} quark mass ratio}'

new_figure = r"""\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{fig_twist_spectrum.pdf}
\caption{A-factor twist spectrum $|\phi(\gamma)|$ (folded to $[0\degr,90\degr]$)
  versus word length for \mfld{m003} (left, green) and \mfld{m006} (right, orange).
  Top panel: $27\degr$--$90\degr$. Bottom panel: $0\degr$--$21\degr$.
  Dashed lines: CKM targets (blue) and PMNS targets (red).
  The spectral gap between $21\degr$ and $27\degr$ is a genuine geometric
  feature of both manifolds.
  Annotated points mark the key coincidences per manifold.}
\label{fig:spectrum}
\end{figure}

"""

if old_anchor in tex:
    tex = tex.replace(old_anchor, new_figure + old_anchor)
    print('Figure block inserted')
else:
    print('WARNING: anchor not found -- searching for alternate')
    # Try the mbmc section label
    alt = r'\label{sec:mbmc}'
    if alt in tex:
        tex = tex.replace(alt, alt + '\n\n' + new_figure)
        print('Figure inserted after sec:mbmc label')
    else:
        print('ERROR: no anchor found -- insert manually')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
