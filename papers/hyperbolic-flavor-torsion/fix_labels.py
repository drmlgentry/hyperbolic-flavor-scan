path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

import re

# 1. Add \label{sec:margulis} to the Margulis subsection
tex = tex.replace(
    r"\subsection{Margulis asymptotics and the exponential decay law}",
    r"\subsection{Margulis asymptotics and the exponential decay law}"
    "\n"
    r"\label{sec:margulis}"
)

# 2. Fix reference to sec:margulis (already uses Section~\ref{sec:margulis} -- now label exists)
# Also fix "this subsection" fallback if present
tex = tex.replace("this subsection", r"Section~\ref{sec:margulis}")

# 3. Ensure Proposition with label is present
# Check if proposition environment exists
if r"\begin{proposition}" not in tex:
    # Add proposition before the Margulis section
    old_margulis_intro = r"\subsection{Margulis asymptotics and the exponential decay law}"
    new_with_prop = (
        r"\subsection{Floors vanish: consequence of the Ambient PGT}"
        "\n"
        r"\label{sec:floors_zero}"
        "\n\n"
        r"\begin{proposition}"
        "\n"
        r"\label{prop:floors_zero}"
        "\n"
        r"Let $M$ be a compact hyperbolic 3-manifold with"
        "\n"
        r"$H_1(M,\mathbb{Z})=\mathbb{Z}/5$."
        "\n"
        r"For each homology class $k\in\mathbb{Z}/5$,"
        "\n"
        r"$\phi_{\rm floor}^{(k)}(L)\to 0$ as $L\to\infty$."
        "\n"
        r"\end{proposition}"
        "\n\n"
        r"\begin{proof}"
        "\n"
        r"By the Ambient Prime Geodesic Theorem~\cite{BalogBiggs23} and the"
        "\n"
        r"equidistribution of homology classes~\cite{Sunada85}, for any"
        "\n"
        r"$\epsilon>0$ the count of class-$k$ geodesics with $\phi_{\rm fold}<\epsilon$"
        "\n"
        r"and length $\leq T$ grows as $(4\epsilon/\pi)\cdot e^{2T}/(10T)\to\infty$."
        "\n"
        r"Hence every class eventually contains geodesics with"
        "\n"
        r"$\phi_{\rm fold}<\epsilon$, so $\phi_{\rm floor}^{(k)}(L)\to 0$."
        "\n"
        r"\end{proof}"
        "\n\n"
        r"\begin{remark}"
        "\n"
        r"Proposition~\ref{prop:floors_zero} shows convergence to zero but"
        "\n"
        r"says nothing about the rate. The class-dependent rates $c_k$"
        "\n"
        r"(Table~\ref{tab:rates}) are the new content of this paper."
        "\n"
        r"\end{remark}"
        "\n\n"
        r"\subsection{Margulis asymptotics and the exponential decay law}"
    )
    tex = tex.replace(old_margulis_intro, new_with_prop, 1)
    print("Proposition inserted.")
else:
    # Just ensure label is present
    if r"\label{prop:floors_zero}" not in tex:
        tex = tex.replace(
            r"\begin{proposition}",
            r"\begin{proposition}" + "\n" + r"\label{prop:floors_zero}"
        )
    print("Proposition already present.")

# 4. Fix remaining "exhibit exponential decay" -> "are well described by..."
tex = tex.replace(
    "exhibit exponential decay",
    "are well described by exponential decay over the observed range"
)

# 5. Fix "stable across lengths >= 10"
tex = tex.replace(
    r"approximately stable for the largest lengths examined ($L=10,11,12$).",
    r"approximately stable for the largest word lengths examined ($L=10,11,12$);"
    "\n"
    r"we make no asymptotic claim beyond the range studied."
)

# 6. Fix "systole of M" refinement
tex = tex.replace(
    "bounded below by the systole $\\mathrm{sys}(\\mfld{m006})>0$",
    "bounded below by the systole of $M$"
)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("All fixes applied.")
