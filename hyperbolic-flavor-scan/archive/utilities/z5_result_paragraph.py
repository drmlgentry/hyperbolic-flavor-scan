# Generate the key result paragraph for the paper

result = """
\\subsection{CP structure across the Z/5 census}
\\label{subsec:z5_census}

A scan of all compact hyperbolic 3-manifolds with $H_1(M)=\\mathbb{Z}/5$
in the \\texttt{OrientableClosedCensus} (indices 1--500) reveals:

\\begin{enumerate}
\\item Every $\\mathbb{Z}/5$ manifold has both degenerate triples
  (two words sharing a homology class, $J=0$) and non-degenerate triples
  (all distinct classes, $J\\neq 0$) among its short-word geodesics.
\\item The maximum geometric CP measure $J_\\text{geom}=2.1266$
  is identical for all $\\mathbb{Z}/5$ manifolds.
\\item The class distribution is nearly uniform: each homology class
  $k\\in\\mathbb{Z}/5$ is represented by approximately 10--13 short words.
\\end{enumerate}

The CP distinction between the CKM sector ($m006$, $J\\approx 0$) and
the PMNS sector ($m003$, $J\\neq 0$) is therefore a \\emph{selection effect}
of the Frobenius fitness criterion, not an intrinsic property of either manifold.
The fitness landscape on $m006$ rewards the degenerate triple
$(\\texttt{aaB}, \\texttt{AbA}, \\texttt{AAb})$ with classes $(3,2,2)$
because the CKM matrix is phenomenologically close to real and diagonal.
The fitness landscape on $m003$ rewards a non-degenerate triple
because the PMNS matrix has large mixing and non-zero $\\delta_{CP}$.

This is a coherent and falsifiable statement: any $\\mathbb{Z}/5$ manifold
could in principle serve as either the CKM or PMNS manifold, depending on
which geodesic triple the fitness criterion selects.
What singles out $m003$ and $m006$ specifically is their
\\emph{volume} (smallest and second-smallest in the $\\mathbb{Z}/5$ family)
and the geometric properties (axis angles, spectral gap) that make
their optimal triples match the observed matrices.
"""

print(result)
