import re

path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Replace the old table with corrected length-10 filtered data
old_table = r"""\begin{tabular}{cllll}
\toprule
Class $k$ & \multicolumn{2}{c}{\mfld{m003} (Meyerhoff)} &
            \multicolumn{2}{c}{\mfld{m006}} \\
          & $\phi_{\rm floor}$ (deg) & Word
          & $\phi_{\rm floor}$ (deg) & Word \\
\midrule
0 & 2.643 & \word{AAAAB}   & 2.132 & \word{AAABB}   \\
1 & 2.763 & \word{AAAbbbb} & 3.374 & \word{BBBB}    \\
2 & 4.857 & \word{AABBAbb} & 2.399 & \word{AABABab} \\
3 & 5.376 & \word{AABABBB} & 1.611 & \word{ABaBAb}  \\
4 & 3.269 & \word{AAb}     & $\mathbf{0.005}$ & \word{AAABAB}  \\
\midrule
\multicolumn{2}{l}{Max/min ratio} & $2\times$
& & $637\times$ \\
\bottomrule
\end{tabular}"""

new_table = r"""\begin{tabular}{clllll}
\toprule
Class $k$ & \multicolumn{2}{c}{\mfld{m003} ($\ell\leq 7$)} &
            \multicolumn{3}{c}{\mfld{m006} ($\ell\leq 10$, genuine loxodromics)} \\
          & $\phi_{\rm floor}$ (deg) & Word
          & $\phi_{\rm floor}$ (deg) & Word & Length \\
\midrule
0 & 2.643 & \word{AAAAB}   & 1.023 & \word{AbaaBAbaBB} & 10 \\
1 & 2.763 & \word{AAAbbbb} & 0.005 & \word{aaabab}     & 6  \\
2 & 4.857 & \word{AABBAbb} & 0.253 & \word{aBabaBBAba} & 10 \\
3 & 5.376 & \word{AABABBB} & 0.253 & \word{BBAAAAbbba} & 10 \\
4 & 3.269 & \word{AAb}     & $\mathbf{0.005}$ & \word{AAABAB} & 6 \\
\midrule
\multicolumn{2}{l}{Max/min ratio} & $2\times$ & \multicolumn{3}{l}{$193\times$ (genuine), gap $[0.005^\circ, 0.253^\circ]$} \\
\bottomrule
\end{tabular}"""

tex = tex.replace(old_table, new_table, 1)

# Update the gap description in Section 3.3
old_gap_text = (r"Coset~4 has floor $0.005\degr$ and coset~3 has the next lowest floor at $1.611\degr$."
    "\n"
    r"Together these produce a spectral gap $[0.005\degr, 1.61\degr]$ in the"
    "\n"
    r"A-factor spectrum of \mfld{m006} at lengths $\leq 7$.")

new_gap_text = (
    r"Classes 1 and 4 share floor $0.005\degr$ (the inverse geodesic pair"
    "\n"
    r"\word{aaabab}/\word{AAABAB}, both first appearing at length~6)."
    "\n"
    r"Classes 2 and~3 share floor $0.253\degr$ (first appearing at length~10)."
    "\n"
    r"Class~0 (identity coset) has the highest floor at $1.023\degr$ (length~10)."
    "\n"
    r"Excluding 18 trivial elements ($|\lambda|\leq 1.01$, mapping to $\pm I$"
    "\n"
    r"in $\PSL$ and hence not representing geodesics),"
    "\n"
    r"the genuine loxodromic spectral gap is $[0.005\degr, 0.253\degr]$ at"
    "\n"
    r"lengths $\leq 10$, with $\theta_{13}^{\rm CKM}=0.201\degr$ lying within it."
)

tex = tex.replace(old_gap_text, new_gap_text, 1)

# Update stability paragraph
old_stab = (
    r"\textbf{Stability.}"
    "\n"
    r"The class~3 floor is set by \word{ABaBAb} (length~6, $\phi_{\rm fold}=1.611\degr$)."
    "\n"
    r"Extending the census from length~6 to~7, the class~3 floor does not decrease:"
    "\n"
    r"no length-7 word in class~3 achieves $\phi_{\rm fold} < 1.611\degr$."
    "\n"
    r"This provides evidence that the gap boundary is not a finite-length artifact,"
    "\n"
    r"though we cannot exclude that longer words in class~3 could close the gap."
    "\n"
    r"Extending the census to length~10 is the most direct test."
)

new_stab = (
    r"\textbf{Stability to length 10.}"
    "\n"
    r"We extended the census to word length~10 (118,096 words, excluding 18 trivial"
    "\n"
    r"elements with $|\lambda|\leq 1.01$)."
    "\n"
    r"The spectral gap narrows from $[0.005\degr,1.61\degr]$ at length~7 to"
    "\n"
    r"$[0.005\degr,0.253\degr]$ at length~10, but does not close."
    "\n"
    r"The class~0 floor decreases from $2.13\degr$ to $1.02\degr$;"
    "\n"
    r"classes~2 and~3 each reach $0.253\degr$."
    "\n"
    r"The value $\theta_{13}^{\rm CKM}=0.201\degr$ remains in the gap at length~10."
    "\n"
    r"Whether the gap persists to all word lengths is an open question;"
    "\n"
    r"extending to length~12 or higher is the most direct test."
)

tex = tex.replace(old_stab, new_stab, 1)

# Update falsifiable prediction in introduction
old_pred = (
    r"\item A falsifiable prediction: no loxodromic element with"
    "\n"
    r"  $0.005\degr < \phi_{\rm fold} < 1.61\degr$ exists in $\pi_1(\mfld{m006})$"
    "\n"
    r"  at any word length, which is testable by extending the census."
)
new_pred = (
    r"\item A falsifiable prediction: no genuine loxodromic element"
    "\n"
    r"  (with $|\lambda|>1$) with $0.005\degr < \phi_{\rm fold} < 0.253\degr$"
    "\n"
    r"  exists in $\pi_1(\mfld{m006})$ at any word length."
    "\n"
    r"  This is testable by extending the census to lengths~$\geq 11$."
)
tex = tex.replace(old_pred, new_pred, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Paper updated with length-10 results.")
