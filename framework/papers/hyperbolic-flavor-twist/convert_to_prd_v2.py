"""
convert_to_prd_v2.py - Convert to revtex4-2, replace longtable with tabular
"""
import re

with open(r"C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex",
          encoding="utf-8") as f:
    content = f.read()

new_preamble = r"""\documentclass[aps,prd,preprint,longbibliography]{revtex4-2}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{booktabs}
\usepackage{array}
\usepackage{microtype}
\usepackage{xcolor}
\usepackage{graphicx}
\usepackage{hyperref}

\newcommand{\PSL}{\mathrm{PSL}(2,\mathbb{C})}
\newcommand{\SLtwo}{\mathrm{SL}(2,\mathbb{C})}
\newcommand{\pioneM}{\pi_1(M)}
\newcommand{\mfld}[1]{\texttt{#1}}
\newcommand{\word}[1]{\texttt{#1}}
\newcommand{\degr}{^{\circ}}
\newcommand{\PDG}{\textsc{pdg}}

\begin{document}

\title{Twist Angle Spectrum of Hyperbolic Holonomy:\\
Encoding of Standard Model Flavor Parameters}

\author{Marvin L.~Gentry}
\affiliation{Independent Researcher, Seattle, WA, USA}
\email{drmlgentry@protonmail.com}
\date{\today}

\begin{abstract}"""

abs_match = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',
                      content, re.DOTALL)
abstract_body = abs_match.group(1).strip() if abs_match else "ABSTRACT MISSING"
new_preamble += "\n" + abstract_body + "\n\\end{abstract}\n\n\\maketitle\n"

body_match = re.search(r'(\\section\{Introduction\}.*?)\\end\{document\}',
                       content, re.DOTALL)
body = body_match.group(1).strip() if body_match else "BODY MISSING"

# Fix acknowledgments
body = re.sub(r'\\section\*\{Acknowledge?ments?\}', r'\\acknowledgments', body)

# Replace longtable with table+tabular
# longtable header: \begin{longtable}{spec}
body = re.sub(r'\\begin\{longtable\}\{[^}]*\}', 
              r'\\begin{table}[htbp]\n\\centering\n\\begin{tabular}{llllll}', body)

# Remove longtable-specific commands
body = re.sub(r'\\endfirsthead\s*', '', body)
body = re.sub(r'\\endhead\s*', '', body)  
body = re.sub(r'\\endfoot\s*', '', body)
body = re.sub(r'\\endlastfoot\s*', '', body)
body = re.sub(r'\\midrule\\multicolumn\{[^}]*\}\{r\}\{\\small continued.*?\}\\\\\s*', '', body)

# Replace \end{longtable} with end of tabular+table
body = re.sub(r'\\end\{longtable\}', r'\\end{tabular}\n\\end{table}', body)

# Move \caption and \label before \begin{tabular} (revtex style)
# Find table blocks and reorder caption/label
def fix_table_caption(m):
    block = m.group(0)
    # Extract caption and label
    cap = re.search(r'(\\caption\{.*?\}(?:\s*\\label\{.*?\})?)', block, re.DOTALL)
    if cap:
        cap_text = cap.group(1)
        block = block.replace(cap_text, '', 1)
        # Insert after \begin{table}[...]
        block = re.sub(r'(\\begin\{table\}[^\n]*\n)', 
                       r'\1' + cap_text + '\n', block)
    return block

body = re.sub(r'\\begin\{table\}.*?\\end\{table\}', 
              fix_table_caption, body, flags=re.DOTALL)

# Fix column specs that revtex dislikes -- replace p{} with l
body = re.sub(r'p\{[^}]+\}', 'l', body)
body = re.sub(r'@\{\}', '', body)
body = re.sub(r'@\{\}', '', body)

# Remove \journal command
body = re.sub(r'\\journal\{[^}]*\}\s*', '', body)

final = new_preamble + "\n\n" + body + "\n\n\\end{document}\n"

out = r"C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-prd.tex"
with open(out, "w", encoding="utf-8") as f:
    f.write(final)
print(f"Written: {out}  ({len(final.splitlines())} lines)")
