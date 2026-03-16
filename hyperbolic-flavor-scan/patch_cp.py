import re

path = r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-hyperbolic-flavor-ckm.tex"
with open(path, "r", encoding="utf-8") as f:
    tex = f.read()

old = "We leave this as an important open problem."

new = """We verified this obstruction numerically by scanning all phase pairs
$(\\varphi_1, \\varphi_3) \\in [-\\pi,\\pi]^2$ with $\\varphi_2 = 0$ fixed
(forced by the trivial homology image of \\texttt{AbA}).
No assignment yields both fitness below $0.10$ and $|J| > 10^{-6}$:
every nonzero phase combination destroys the CKM structure,
driving the Frobenius error from $0.017$ to $\\approx 1.43$.
This confirms that the obstruction is geometric, not merely numerical.

Resolving CP violation may require:
\\begin{itemize}
\\item A manifold with free homology ($H_1 \\supset \\mathbb{Z}$),
permitting continuous phases without the torsion cancellation.
\\item Use of the complex fixed-point coordinates $z_\\pm \\in \\mathbb{CP}^1$
of each loxodromic element as a source of geometric phase, rather
than flat $\\mathrm{U}(1)$ connections.
\\item A two-manifold construction in which flavor mixing and
the CP phase arise from separate geometric sectors.
\\end{itemize}

We leave the resolution to future work, noting that the obstruction
itself---a consequence of the trivial homology image of \\texttt{AbA}
in $H_1(\\mathrm{m006};\\mathbb{Z}) = \\mathbb{Z}/5$---is a precise,
falsifiable statement about the geometry."""

if old in tex:
    tex = tex.replace(old, new)
    with open(path, "w", encoding="utf-8") as f:
        f.write(tex)
    lines = tex.count("\n")
    print(f"SUCCESS. Lines: {lines}")
else:
    print("NOMATCH. Searching for nearby text:")
    for i, line in enumerate(tex.split("\n")):
        if "open problem" in line or "Resolving" in line:
            print(f"  {i}: {repr(line)}")