"""
convert_to_prd.py - Convert NPB elsarticle tex to revtex4-2/prd format
"""
import re

with open(r"C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex",
          encoding="utf-8") as f:
    content = f.read()

# ── New preamble (replaces everything before \begin{document}) ────
new_preamble = r"""\documentclass[aps,prd,preprint,longbibliography]{revtex4-2}

\usepackage{amsmath,amssymb,amsfonts}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{graphicx}
\usepackage{hyperref}

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

\title{Twist Angle Spectrum of Hyperbolic Holonomy:\\
Encoding of Standard Model Flavor Parameters}

\author{Marvin L.~Gentry}
\affiliation{Independent Researcher, Seattle, WA, USA}
\email{drmlgentry@protonmail.com}

\date{\today}

\begin{abstract}"""

# Find the abstract content from the NPB version
abs_match = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',
                      content, re.DOTALL)
abstract_body = abs_match.group(1).strip() if abs_match else "ABSTRACT MISSING"

new_preamble += "\n" + abstract_body + "\n\\end{abstract}\n\n\\maketitle\n"

# Extract everything from first \section to \end{document}
body_match = re.search(r'(\\section\{Introduction\}.*?)\\end\{document\}',
                       content, re.DOTALL)
if not body_match:
    # Try from any section
    body_match = re.search(r'(\\section.*?)\\end\{document\}', content, re.DOTALL)

body = body_match.group(1).strip() if body_match else "BODY MISSING"

# Fix acknowledgments: elsarticle uses \section*{} revtex uses \acknowledgments
body = body.replace(r'\section*{Acknowledgements}', r'\acknowledgments')
body = body.replace(r'\section*{Acknowledgments}', r'\acknowledgments')

# Remove \journal{} command if present
body = re.sub(r'\\journal\{[^}]*\}\s*', '', body)

# Assemble final document
final = new_preamble + "\n\n" + body + "\n\n\\end{document}\n"

out = r"C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-prd.tex"
with open(out, "w", encoding="utf-8") as f:
    f.write(final)

print(f"Written: {out}")
print(f"Lines: {len(final.splitlines())}")
