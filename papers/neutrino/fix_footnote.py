path = "gentry-neutrino.tex"
with open(path, encoding="utf-8") as f:
    src = f.read()

old = r"""\end{tabular}
\footnotetext{The other bases produce different integer triples
$(q_1,q_2,q_3)$ than the $\phi$ triple; they are not equivalent
to the $\phi$ solution.}"""

new = r"""\end{tabular}
\smallskip\par
{\small \textit{Note:} The other bases produce different integer
triples $(q_1,q_2,q_3)$ than the $\phi$ triple; they are not
equivalent to the $\phi$ solution.}"""

if old in src:
    src = src.replace(old, new)
    print("Fixed footnote")
else:
    print("Pattern not found - check manually")
    # Show context around footnotetext
    idx = src.find("footnotetext")
    if idx >= 0:
        print(repr(src[idx-50:idx+150]))

with open(path, "w", encoding="utf-8") as f:
    f.write(src)
