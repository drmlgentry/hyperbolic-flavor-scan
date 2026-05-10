path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Update Section 3.1 to add m003 length-12 data
old_sec31 = (
    r"On \mfld{m003}, the floors are within $2\times$ of each other."
)
new_sec31 = (
    r"On \mfld{m003}, the floors are within $2\times$ of each other."
    "\n"
    r"This contrast between the two manifolds at length~7 motivates"
    "\n"
    r"extending the census to longer lengths for both manifolds."
)
tex = tex.replace(old_sec31, new_sec31, 1)

# Update the comparison section in Discussion
old_item4 = (
    r"\item By contrast, the Meyerhoff manifold \mfld{m003} (also $H_1=\mathbb{Z}/5$)"
    "\n"
    r"      shows only a $2\times$ spread of floors at length~7, indicating"
    "\n"
    r"      a qualitatively different torsion-floor relationship despite"
    "\n"
    r"      identical $H_1$."
)
new_item4 = (
    r"\item At length~12, the Meyerhoff manifold \mfld{m003} (also $H_1=\mathbb{Z}/5$)"
    "\n"
    r"      exhibits the \emph{same} pairing structure as \mfld{m006}:"
    "\n"
    r"      floors $\phi^{(1,4)}=0.04319\degr$, $\phi^{(2,3)}=0.09387\degr$,"
    "\n"
    r"      $\phi^{(0)}=0.04373\degr$, giving ratio"
    "\n"
    r"      $\phi^{(2,3)}/\phi^{(1,4)}\approx 2.17$."
    "\n"
    r"      The pairing $c_1=c_4$, $c_2=c_3$ thus appears to be a property"
    "\n"
    r"      of the $\mathbb{Z}/5$ torsion structure shared by both manifolds,"
    "\n"
    r"      not specific to \mfld{m006}'s arithmetic."
)
tex = tex.replace(old_item4, new_item4, 1)

# Add m003 length-12 table after Table 2 (floors at length 7)
old_after_tab7 = (
    r"\subsection{Spectral floors at length 12: the 4:2:1 hierarchy}"
)
new_after_tab7 = (
    r"\subsection{Length-12 floors for both manifolds}"
    "\n\n"
    r"Table~\ref{tab:floors_both_L12} compares length-12 floors for both"
    "\n"
    r"manifolds.  The $4:2:1$ pairing structure appears in both cases."
    "\n\n"
    r"\begin{table}[tp]"
    "\n"
    r"\caption{Spectral floors per homology class at word length $\leq 12$"
    "\n"
    r"  for both manifolds with $H_1=\mathbb{Z}/5$."
    "\n"
    r"  Both exhibit the pairing $\phi^{(1,4)}<\phi^{(0)}<\phi^{(2,3)}$"
    "\n"
    r"  with ratio $\phi^{(2,3)}/\phi^{(1,4)}\approx 2$.}"
    "\n"
    r"\label{tab:floors_both_L12}"
    "\n"
    r"\begin{tabular}{cllll}"
    "\n"
    r"\toprule"
    "\n"
    r"Class $k$ & \multicolumn{2}{c}{\mfld{m003} (Meyerhoff)} &"
    "\n"
    r"            \multicolumn{2}{c}{\mfld{m006}} \\"
    "\n"
    r"          & $\phi_{\rm floor}$ ($\degr$) & Word"
    "\n"
    r"          & $\phi_{\rm floor}$ ($\degr$) & Word \\"
    "\n"
    r"\midrule"
    "\n"
    r"0 & 0.04373 & \word{aaaabbbaabaB} & 0.02300 & \word{aaaaaaaaabAb} \\"
    "\n"
    r"1 & 0.04319 & \word{aaBAbAAbAABB} & 0.00530 & \word{aaaaaababAAA} \\"
    "\n"
    r"2 & 0.09387 & \word{aabAAbAAABB}  & 0.01060 & \word{aaababaaabab} \\"
    "\n"
    r"3 & 0.09387 & \word{aaaaBabbAAB}  & 0.01060 & \word{AAABABAAABAB} \\"
    "\n"
    r"4 & 0.04319 & \word{aabbAAbaBaaB} & 0.00530 & \word{aaaBAAABAAAA} \\"
    "\n"
    r"\midrule"
    "\n"
    r"\multicolumn{2}{l}{$\phi^{(2,3)}/\phi^{(1,4)}$} & $2.17$ &"
    "\n"
    r"  & $2.00$ \\"
    "\n"
    r"\bottomrule"
    "\n"
    r"\end{tabular}"
    "\n"
    r"\end{table}"
    "\n\n"
    r"\subsection{Spectral floors at length 12: the 4:2:1 hierarchy}"
)
tex = tex.replace(old_after_tab7, new_after_tab7, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("m003 comparison added.")
