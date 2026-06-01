import re

fig_block = r"""
\begin{figure*}[tp]
\includegraphics[width=\textwidth]{fig_coset_spectrum.pdf}
\caption{A-factor twist spectrum $|\phi(\gamma)|$ (folded to $[0\degr,90\degr]$)
  vs.\ word length for \mfld{m006}, coloured by homology class $k\in\mathbb{Z}/5$.
  The isolated class-4 point (\word{AAABAB}, $\phi_{\rm fold}=0.0053\degr$, length~6)
  is annotated; the grey band marks the spectral gap $[0.005\degr,1.61\degr]$;
  the dashed line marks $\theta_{13}^{\rm CKM}=0.201\degr$.}
\label{fig:spectrum}
\end{figure*}

\begin{figure*}[tp]
\includegraphics[width=\textwidth]{fig_coset_floors.pdf}
\caption{Spectral floors $\phi_{\rm floor}(k)$ per homology class $k\in\mathbb{Z}/5$
  for \mfld{m003} (left, Meyerhoff) and \mfld{m006} (right), from the
  word-length-$\leq 7$ census.
  The $637\times$ asymmetry on \mfld{m006} (class-4 floor: $0.005\degr$;
  next lowest: $1.61\degr$) produces the spectral gap.
  On \mfld{m003} all five floors lie within a factor of $2\times$.}
\label{fig:floors}
\end{figure*}

"""

path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

tex = tex.replace(r"\section{Discussion}", fig_block + r"\section{Discussion}", 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Figures inserted.")
