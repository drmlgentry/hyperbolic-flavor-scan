path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

fig1 = """
\\begin{figure*}[!t]
\\centering
\\includegraphics[width=0.78\\textwidth]{fig_cp_pareto}
\\caption{Pareto frontier in the $({\\mathcal F}, |J|)$ plane for the
complex Borel construction on $m003$ with word triple
\\texttt{Ba/baB/BBa}, as the imaginary weight $\\alpha$ varies
from $0$ to $1.125$. The real construction ($\\alpha = 0$, green dashed line)
achieves $\\mathcal{F} = 0.019$ with $J = 0$. Increasing $\\alpha$
continuously trades fitness for CP phase: at $\\alpha = 0.9$ the fitness
remains $0.020$ while $|J| = 0.013$; at $\\alpha = 1.125$ the Jarlskog
invariant matches the PDG value $|J_{\\rm PDG}| = 0.00966$ (red dashed line)
exactly. Color encodes $\\alpha$.}
\\label{fig:pareto}
\\end{figure*}
"""

fig2 = """
\\begin{figure*}[!t]
\\centering
\\includegraphics[width=0.92\\textwidth]{fig_cp_afactor}
\\caption{A-factor twist angles $\\phi(\\gamma)$ for the optimal word triples
on $m003$ (PMNS, left) and $m006$ (CKM, right). On $m003$, the three words
\\texttt{aa}, \\texttt{ab}, \\texttt{aB} have three distinct twist angles
($-168.6^\\circ$, $83.7^\\circ$, $95.7^\\circ$), producing a non-degenerate
$A$-factor and $J_{\\rm PMNS} \\approx -0.010$.
On $m006$, all three words \\texttt{aaB}, \\texttt{AbA}, \\texttt{AAb}
share the same twist angle $\\phi \\approx 92.5^\\circ \\approx \\pi/2$,
producing a degenerate $A$-factor and the observed suppression
$J_{\\rm CKM} \\approx 3\\times10^{-5} \\ll J_{\\rm PMNS}$.
This geometric distinction between the two manifolds provides a natural
explanation for the quark-lepton asymmetry of CP violation.}
\\label{fig:afactor}
\\end{figure*}
"""

# Fig 1 goes before the complex results subsection
marker1 = r'\subsection{Numerical results}'
if marker1 in tex:
    tex = tex.replace(marker1, fig1 + '\n' + marker1, 1)
    print('Fig 1 inserted before numerical results')
else:
    print('Fig 1 marker not found')

# Fig 2 goes before the quark-lepton asymmetry subsection
marker2 = r'\subsection{Quark-lepton asymmetry of CP violation'
if marker2 in tex:
    tex = tex.replace(marker2, fig2 + '\n' + marker2, 1)
    print('Fig 2 inserted before quark-lepton subsection')
else:
    print('Fig 2 marker not found')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
