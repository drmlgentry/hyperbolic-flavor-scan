import re

# Read the current file
path = r"C:\dev\framework\papers\holonomy-cp\gentry-holonomy-cp.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Fix stale IDs
tex = tex.replace("es2026mar11_966", "DQ14014")
tex = tex.replace("es2026mar13_942", "DQ14050")

# Fix abstract
tex = tex.replace(
    "and illustrate the mechanism using the Weeks manifold.",
    r"and illustrate the mechanism using the Meyerhoff manifold ($\mathrm{OrientableClosedCensus}[1]$, $\mathrm{vol}=0.9814$, $H_1=\mathbb{Z}/5$)."
)

# Fix section title
tex = tex.replace(
    r"\section{Example: The Weeks Manifold}",
    r"\section{Example: The Meyerhoff Manifold}"
)

# Fix example section body
old_example = r"""The Weeks manifold $M_W$ is the closed hyperbolic $3$-manifold
of smallest known volume.  Its first homology group satisfies
\[
H_1(M_W;\mathbb{Z}) \cong (\mathbb{Z}/5)^3.
\]
Hence there exist nontrivial characters
\[
\chi : H_1(M_W;\mathbb{Z}) \to U(1)
\]
sending generators to fifth roots of unity.
Choose $\rho$ corresponding to such a character and let $\gamma$ be a
geodesic representing a generator of a $\mathbb{Z}/5$ factor.
One can check that $M_W$ admits an orientation-reversing isometry
that preserves the conjugacy class of $\gamma$; for instance, the
involution described in~\cite{GoodmanEtAl} acts by inversion on
homology, mapping $\chi$ to its complex conjugate.
Therefore
\[
W_\rho(\gamma) = e^{2\pi i/5},
\]
which is $\neq \pm 1$ and provides a concrete geometric CP-odd phase."""

new_example = r"""The Meyerhoff manifold $M = \mathrm{OrientableClosedCensus}[1]$
(SnapPy notation~\cite{SnapPy}, $\mathrm{vol}=0.9814$) is the compact
oriented hyperbolic $3$-manifold of second-smallest known volume.
It is proven arithmetic with invariant trace field of discriminant
$-283$~\cite{Chinburg98}.  Its first homology group satisfies
\[
H_1(M;\mathbb{Z}) \cong \mathbb{Z}/5.
\]
Hence there exist nontrivial characters
\[
\chi : H_1(M;\mathbb{Z}) \to U(1)
\]
sending the generator to a fifth root of unity.
Choose $\rho$ corresponding to such a character and let $\gamma$ be a
geodesic representing the generator of $\mathbb{Z}/5$.
The manifold $M$ admits an orientation-reversing isometry acting by
inversion on $H_1(M;\mathbb{Z})=\mathbb{Z}/5$, mapping $\chi$ to its
complex conjugate; this is consistent with complex conjugation in the
invariant trace field acting on the holonomy representation~\cite{MaclachlanReid}.
Therefore
\[
W_\rho(\gamma) = e^{2\pi i/5},
\]
which is $\neq \pm 1$ and provides a concrete geometric CP-odd phase."""

tex = tex.replace(old_example, new_example, 1)

# Fix conclusion
tex = tex.replace(
    "The Weeks manifold provides an explicit example realising a fifth root\nof unity as a CP-odd phase.",
    "The Meyerhoff manifold provides an explicit example realising a fifth root\nof unity as a CP-odd phase."
)

# Fix bibliography: replace GoodmanEtAl with Chinburg98 and MaclachlanReid
tex = tex.replace(
    r"""\bibitem{GoodmanEtAl}
O.~Goodman, D.~Heard, and C.~Hodgson,
``Commensurators of cusped hyperbolic manifolds,''
\emph{Exp. Math.} \textbf{17} (2008), 283.""",
    r"""\bibitem{Chinburg98}
T.~Chinburg, E.~Friedman, K.~N.~Jones, and A.~W.~Reid,
``The arithmetic hyperbolic 3-manifold of smallest volume,''
\emph{Ann. Scuola Norm. Sup. Pisa} \textbf{30}, 1--40 (2001).

\bibitem{MaclachlanReid}
C.~Maclachlan and A.~W.~Reid,
\emph{The Arithmetic of Hyperbolic 3-Manifolds},
Springer, New York (2003)."""
)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)

# Verify
checks = ["Meyerhoff", "Weeks manifold", "DQ14014", "DQ14050", "Chinburg98", "MaclachlanReid"]
for c in checks:
    count = tex.count(c)
    print(f"  {c}: {count}")
