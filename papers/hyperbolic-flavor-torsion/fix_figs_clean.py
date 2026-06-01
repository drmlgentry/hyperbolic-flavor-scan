path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

import re

# Remove ALL existing figure* blocks -- we will reinsert clean ones
tex = re.sub(r'\\begin\{figure\*\}.*?\\end\{figure\*\}', '', tex, flags=re.DOTALL)

# Clean figure block to insert before \section{Discussion}
figs = (
    "\n\\begin{figure*}[tp]\n"
    "\\includegraphics[width=\\textwidth]{fig_coset_spectrum.pdf}\n"
    "\\caption{A-factor twist spectrum $|\\phi(\\gamma)|$ vs.\\ word length "
    "for \\mfld{m006}, coloured by homology class $k\\in\\mathbb{Z}/5$. "
    "The isolated class-4 point (\\word{AAABAB}, "
    "$\\phi_{\\rm fold}=0.0053^{\\circ}$, length~6) is annotated; "
    "the grey band marks the spectral gap $[0.005^{\\circ},1.61^{\\circ}]$; "
    "the dashed line marks $\\theta_{13}^{\\rm CKM}=0.201^{\\circ}$.}\n"
    "\\label{fig:spectrum}\n"
    "\\end{figure*}\n\n"
    "\\begin{figure*}[tp]\n"
    "\\includegraphics[width=\\textwidth]{fig_coset_floors.pdf}\n"
    "\\caption{Spectral floors $\\phi_{\\rm floor}(k)$ per homology class "
    "$k\\in\\mathbb{Z}/5$ for \\mfld{m003} (left, Meyerhoff) and \\mfld{m006} "
    "(right), from the word-length-$\\leq 7$ census. "
    "The $637\\times$ asymmetry on \\mfld{m006} "
    "(class-4 floor: $0.005^{\\circ}$; next: $1.61^{\\circ}$) "
    "produces the spectral gap. "
    "On \\mfld{m003} all floors lie within $2\\times$ of each other.}\n"
    "\\label{fig:floors}\n"
    "\\end{figure*}\n\n"
)

tex = tex.replace(r"\section{Discussion}", figs + r"\section{Discussion}", 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Clean figures inserted.")
