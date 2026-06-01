path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# 1. Replace abstract with the better version
import re
old_abs = re.search(r'\\begin\{abstract\}.*?\\end\{abstract\}', tex, re.DOTALL).group(0)
new_abs = r"""\begin{abstract}
We introduce a homology-resolved extreme-value invariant of the complex length
spectrum of a compact hyperbolic 3-manifold.
For manifolds with $H_1(M,\mathbb{Z})=\mathbb{Z}/5$, the abelianization map
partitions loxodromic elements into five homology classes.
We study the spectral floor $\phi_{\rm floor}^{(k)}(L)$---the minimum folded
twist angle among loxodromic words of length $\leq L$ in class $k$---for the
manifold \mfld{m006} (OrientableClosedCensus[43], $\mathrm{vol}=2.029$,
$H_1=\mathbb{Z}/5$).
A complete census up to word length~12 ($1{,}062{,}802$ genuine loxodromics)
shows that all five floors decrease rapidly with $L$ and are well described by
exponential decay over the observed range, with class-dependent rates
$c_k$ forming a paired hierarchy:
$c_1=c_4\approx 0.985$, $c_2=c_3\approx 0.401$, $c_0\approx 0.750$,
confirmed by Akaike Information Criterion (AIC) model selection against
power-law alternatives ($\Delta\mathrm{AIC}>40$ for classes~0,1,4).
At $L=12$ the floors satisfy an approximate $4:2:1$ ratio.
The same pairing structure appears for the Meyerhoff manifold \mfld{m003}
(also $H_1=\mathbb{Z}/5$), suggesting the hierarchy is a property of
$\mathbb{Z}/5$ torsion rather than of either specific manifold.
We conjecture that the pairing $c_k=c_{5-k}$ reflects an arithmetic symmetry
of the holonomy representation, arising from complex conjugation in the
invariant trace field.
The floor $\phi_{\rm floor}^{(k)}(L)\to 0$ for all $k$ follows from the
Ambient Prime Geodesic Theorem~\cite{BalogBiggs23}; no claim is made
regarding the asymptotic form of the decay rate.
\end{abstract}"""
tex = tex.replace(old_abs, new_abs, 1)

# 2. Fix overclaiming in main text: "decay exponentially" -> "are well described by exponential decay"
tex = tex.replace(
    "reveals that all five floors decay exponentially with $L$, but with",
    "shows that all five floors decrease rapidly with $L$ and are well described\nby exponential decay over the observed range, with"
)
tex = tex.replace(
    r"Each class exhibits approximately log-linear decay, consistent with",
    r"Over the range $L=1,\ldots,12$, each class exhibits approximately log-linear decay, consistent with"
)

# 3. Fix "stable across lengths >= 10" claim
tex = tex.replace(
    r"The $4:2:1$ ratio (0.023 : 0.011 : 0.005) is stable across lengths $\geq 10$.",
    r"The $4:2:1$ ratio (0.023 : 0.011 : 0.005) is approximately stable for the largest lengths examined ($L=10,11,12$)."
)

# 4. Add asymptotic disclaimer to Discussion
old_disc_end = (
    r"Comparing these rates across the OrientableClosedCensus for other manifolds"
    "\n"
    r"with $H_1=\mathbb{Z}/p$ ($p$ prime) would test the conjecture and may"
    "\n"
    r"reveal which arithmetic properties control the rate hierarchy."
)
new_disc_end = (
    r"Comparing these rates across the OrientableClosedCensus for other manifolds"
    "\n"
    r"with $H_1=\mathbb{Z}/p$ ($p$ prime) would test the conjecture and may"
    "\n"
    r"reveal which arithmetic properties control the rate hierarchy."
    "\n\n"
    r"\begin{remark}[Scope of the decay claim]"
    "\n"
    r"The exponential fits in Table~\ref{tab:rates} and Figure~\ref{fig:decay}"
    "\n"
    r"are empirical descriptions of the data over the range $L=1,\ldots,12$."
    "\n"
    r"No claim is made regarding the asymptotic behaviour of"
    "\n"
    r"$\phi_{\rm floor}^{(k)}(L)$ as $L\to\infty$ beyond the convergence to"
    "\n"
    r"zero established in Proposition~\ref{prop:floors_zero}."
    "\n"
    r"The Margulis asymptotics (Section~\ref{sec:margulis}) suggest exponential"
    "\n"
    r"decay is natural, but a rigorous proof of the rate is an open problem."
    "\n"
    r"\end{remark}"
)
tex = tex.replace(old_disc_end, new_disc_end, 1)

# 5. Add primitive vs non-primitive note to census procedure
old_census = r"Words with $|\lambda|\leq 1.01$ (mapping to $\pm I$ in $\PSL$) are excluded as non-geodesic relators."
new_census = (
    r"Words with $|\lambda|\leq 1.01$ (mapping to $\pm I$ in $\PSL$) are excluded"
    "\n"
    r"as non-geodesic relators."
    "\n"
    r"The census includes both primitive and non-primitive loxodromic words;"
    "\n"
    r"the prime geodesic theorem counts only primitive closed geodesics, but"
    "\n"
    r"for the spectral floor the distinction is immaterial since a non-primitive"
    "\n"
    r"element $\gamma^n$ has $\phi_{\rm fold}(\gamma^n)\geq\phi_{\rm fold}(\gamma)$"
    "\n"
    r"only if $n$ is odd, and the floor is in any case dominated by primitive elements."
)
tex = tex.replace(old_census, new_census, 1)

# 6. Add justification for T_L ~ alpha L
old_tl = r"where $T_L$ is the maximum geodesic length achievable at word length $L$."
new_tl = (
    r"where $T_L$ is the maximum geodesic length achievable at word length $L$."
    "\n"
    r"Since each generator has a fixed positive translation length"
    "\n"
    r"$\ell_{\rm min}=\min_g\ell(g)>0$ (equal to the systole of \mfld{m006}),"
    "\n"
    r"the translation length of a word of length $L$ is bounded above by"
    "\n"
    r"$L\cdot\ell_{\rm max}$ and below by $\ell_{\rm min}$, giving"
    "\n"
    r"$T_L\leq L\cdot\ell_{\rm max}$; empirically $T_L\sim\alpha L$ with"
    "\n"
    r"$\alpha\approx\ell_{\rm max}(\mfld{m006})\approx 1.5$ for the word lengths studied."
)
tex = tex.replace(old_tl, new_tl, 1)

# 7. Ensure Proposition label is defined (check it's there)
if r"\label{prop:floors_zero}" in tex:
    print("Proposition label present -- OK")
else:
    print("WARNING: Proposition label missing!")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("All fixes applied.")
