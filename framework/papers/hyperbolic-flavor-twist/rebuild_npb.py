import re

path_in  = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist.tex'
path_out = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'

with open(path_in, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Extract abstract text
m = re.search(r'\\abstract\{(.*?)\n\}', tex, re.DOTALL)
abstract = m.group(1).strip() if m else 'ABSTRACT NOT FOUND'

# Extract body: everything from \section{Introduction} to \end{document}
m2 = re.search(r'(\\section\{Introduction\}.*?)\\end\{document\}', tex, re.DOTALL)
body = m2.group(1).strip() if m2 else 'BODY NOT FOUND'

# Fix acknowledgments
body = body.replace(r'\acknowledgments', r'\section*{Acknowledgements}')

out = r"""\documentclass[preprint,12pt,number]{elsarticle}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{graphicx}
\usepackage{lineno}

\newcommand{\PSL}{\mathrm{PSL}(2,\mathbb{C})}
\newcommand{\SLtwo}{\mathrm{SL}(2,\mathbb{C})}
\newcommand{\pioneM}{\pi_1(M)}
\newcommand{\phiA}[1]{\phi(#1)}
\newcommand{\mfld}[1]{\texttt{#1}}
\newcommand{\word}[1]{\texttt{#1}}
\newcommand{\degr}{^{\circ}}
\newcommand{\PDG}{\textsc{pdg}}
\newcolumntype{C}{>{$}c<{$}}

\begin{document}

\journal{Nuclear Physics B}

\begin{frontmatter}

\title{Twist Angle Spectrum of Hyperbolic Holonomy:
Encoding of Standard Model Flavor Parameters}

\author[a]{Marvin L.~Gentry}
\ead{drmlgentry@protonmail.com}
\address[a]{Independent Researcher, Seattle, WA, USA}

\begin{abstract}
""" + abstract + r"""
\end{abstract}

\begin{keyword}
hyperbolic 3-manifolds \sep holonomy \sep loxodromic elements \sep
Standard Model flavor \sep CKM \sep PMNS \sep quark masses \sep
Selberg zeta function \sep Iwasawa decomposition
\end{keyword}

\end{frontmatter}

\linenumbers

""" + body + r"""

\end{document}
"""

with open(path_out, 'w', encoding='utf-8') as f: f.write(out)
print('Clean NPB file written.')
