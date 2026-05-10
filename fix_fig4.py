with open(r"C:\dev\golden_unification\papers\CORE_MASTER_v2.tex", encoding="utf-8") as f:
    tex = f.read()

TIKZ_FIG = r"""
\begin{tikzpicture}[
    node distance=2.0cm,
    box/.style={draw, rounded corners, minimum width=4cm, minimum height=1.2cm, align=center, font=\small},
    topbox/.style={box, fill=yellow!20},
    bluebox/.style={box, fill=blue!10},
    greenbox/.style={box, fill=green!15},
    redbox/.style={box, fill=red!15}
]
\node[topbox] (arith) at (0,4) {Arithmetic structure: $\mathbb{Q}(\sqrt{5})$ unit lattice, $H_1(M)=\mathbb{Z}/5$};
\node[bluebox, below left=of arith] (symm) {Symmetric space \\ $\mathrm{SL}(2,\mathbb{R})/\mathrm{SO}(2)$};
\node[bluebox, below right=of arith] (hyp) {Hyperbolic 3-manifolds \\ $m003,\, m006$};
\node[greenbox, below=of symm] (mass) {Mass hierarchy \\ $m \sim m_e \phi^{q/4}$};
\node[bluebox, below=of arith, yshift=-1.5cm] (cp) {CP chirality theorem \\[3pt] Amphicheiral $\Leftrightarrow \exists\, g$ with $g_*=-1$ on $H_1$ \\ Chiral $\Leftrightarrow$ no such $g$};
\node[greenbox, below=of hyp] (iwasawa) {Iwasawa decomposition \\ $\mathrm{PSL}(2,\mathbb{C}) = KAN$};
\node[redbox, below=of mass] (ckm) {CKM mixing \\ ($K$-factor, $m006$, $p<10^{-4}$)};
\node[redbox, below=of cp] (phase) {CP phase \\ $e^{2\pi ik/5}$, $\delta_{CP}=203.5^\circ$};
\node[redbox, below=of iwasawa] (pmns) {PMNS mixing \\ ($N$-factor, $m003$, $p=0.0002$)};
\draw[->] (arith) -- (symm);
\draw[->] (arith) -- (hyp);
\draw[->] (symm) -- (mass);
\draw[->] (hyp) -- (iwasawa);
\draw[->] (mass) -- (ckm);
\draw[->] (cp) -- (phase);
\draw[->] (iwasawa) -- (pmns);
\draw[->] (symm) -- (cp);
\draw[->] (hyp) -- (cp);
\end{tikzpicture}
"""

# Replace the includegraphics line for fig4
old = r"\includegraphics[width=0.85\textwidth]{fig4_framework_diagram}"
tex = tex.replace(old, TIKZ_FIG, 1)

with open(r"C:\dev\golden_unification\papers\CORE_MASTER_v2.tex", "w", encoding="utf-8") as f:
    f.write(tex)
print("TikZ inserted:", TIKZ_FIG[:50] in tex)
print("Lines:", len(tex.splitlines()))
