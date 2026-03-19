path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# ── 1. Switch to CMP/Springer format ──────────────────────────────
old_preamble = r"""\documentclass[aps,prd,twocolumn,longbibliography]{revtex4-2}
\usepackage{amsmath,amssymb,amsfonts,booktabs,array,microtype,xcolor,graphicx,hyperref}

\newcommand{\PSL}{\mathrm{PSL}(2,\mathbb{C})}
\newcommand{\mfld}[1]{\texttt{#1}}
\newcommand{\word}[1]{\texttt{#1}}
\newcommand{\degr}{^{\circ}}
\newcommand{\ZZ}{\mathbb{Z}}

\begin{document}

\title{Homology Class Asymmetry and Spectral Floors in the A-Factor Spectrum\\
of Compact Hyperbolic 3-Manifolds}

\author{Marvin L.~Gentry}
\affiliation{Independent Researcher, Seattle, WA, USA}
\email{drmlgentry@protonmail.com}
\date{\today}"""

new_preamble = r"""\documentclass[smallextended]{svjour3}
\usepackage{amsmath,amssymb,amsfonts,booktabs,array,microtype,xcolor,graphicx,hyperref}
\smartqed

\newcommand{\PSL}{\mathrm{PSL}(2,\mathbb{C})}
\newcommand{\mfld}[1]{\texttt{#1}}
\newcommand{\word}[1]{\texttt{#1}}
\newcommand{\degr}{^{\circ}}
\newcommand{\ZZ}{\mathbb{Z}}

\journalname{Communications in Mathematical Physics}

\begin{document}

\title{Homology Class Asymmetry and Exponential Decay of Spectral Floors
in the Loxodromic Twist Spectrum of Compact Hyperbolic 3-Manifolds}

\author{Marvin L.~Gentry}
\institute{Independent Researcher, Seattle, WA, USA \\
\email{drmlgentry@protonmail.com}}
\date{\today}"""

tex = tex.replace(old_preamble, new_preamble, 1)

# ── 2. Fix maketitle block ─────────────────────────────────────────
tex = tex.replace(r"\maketitle", r"\maketitle", 1)  # already fine

# ── 3. Fix acknowledgments for Springer ───────────────────────────
tex = tex.replace(r"\acknowledgments", r"\begin{acknowledgements}")
tex = tex.replace(
    r"Computations used \textsc{SnapPy}~\cite{SnapPy}, NumPy, and Pandas.",
    r"Computations used \textsc{SnapPy}~\cite{SnapPy}, NumPy, and Pandas."
    "\n" + r"\end{acknowledgements}"
)

# ── 4. Add conjecture to Discussion ───────────────────────────────
old_disc_end = (
    r"\item The exponential decay rates $c_k$ extracted from the log-linear plot"
    "\n"
    r"      provide quantitative data that can be compared with predictions from"
    "\n"
    r"      the length spectrum and the Selberg zeta function.  A full analysis"
    "\n"
    r"      of the asymptotic distribution of twist angles is left to future work."
    "\n"
    r"\end{enumerate}"
)

new_disc_end = (
    r"\item The exponential decay rates $c_k$ extracted from the log-linear plot"
    "\n"
    r"      provide quantitative data that can be compared with predictions from"
    "\n"
    r"      the length spectrum and the Selberg zeta function.  A full analysis"
    "\n"
    r"      of the asymptotic distribution of twist angles is left to future work."
    "\n"
    r"\end{enumerate}"
    "\n\n"
    r"\begin{conjecture}"
    "\n"
    r"The decay rates $c_k$ are determined by the arithmetic structure of"
    "\n"
    r"$\pi_1(\mfld{m006})$. The pairing $c_1=c_4$ and $c_2=c_3$ reflects an"
    "\n"
    r"outer involution $\rho\mapsto\bar\rho$ of the holonomy representation"
    "\n"
    r"$\rho:\pi_1(\mfld{m006})\to\PSL$, which exchanges classes $k\leftrightarrow 5-k$"
    "\n"
    r"(mod~5). The unpaired class $k=0$ has a distinct rate $c_0\approx 0.75$."
    "\n"
    r"More precisely, manifolds with the same invariant trace field and"
    "\n"
    r"quaternion algebra as \mfld{m006} have the same decay rate hierarchy"
    "\n"
    r"$(c_0,c_1,c_2,c_3,c_4)$."
    "\n"
    r"\end{conjecture}"
    "\n\n"
    r"A proof would require computing the invariant trace field of \mfld{m006}"
    "\n"
    r"via the \textsc{Snap} package~\cite{Snap} and relating the arithmetic"
    "\n"
    r"involution to the observed class pairing."
)

tex = tex.replace(old_disc_end, new_disc_end, 1)

# ── 5. Add decay rate table ────────────────────────────────────────
old_decay_eq = (
    r"with slopes $c_k$ that differ by factors of roughly $2$ between tiers."
    "\n"
    r"This exponential approach to zero is expected from generic hyperbolic dynamics"
    "\n"
    r"(closed geodesics become dense), but the class-dependent rates encode"
    "\n"
    r"the torsion structure."
)

new_decay_eq = (
    r"with slopes $c_k$ that differ by factors of roughly $2$ between tiers."
    "\n"
    r"Table~\ref{tab:rates} gives the fitted rates and AIC scores confirming"
    "\n"
    r"exponential over power-law decay for all five classes."
    "\n\n"
    r"\begin{table}[tp]"
    "\n"
    r"\caption{Exponential decay rates $c_k$ from log-linear fits to"
    "\n"
    r"$\phi_{\rm floor}^{(k)}(L)\sim e^{-c_k L}$, with AIC comparison"
    "\n"
    r"against power-law fits. Lower AIC = better fit.}"
    "\n"
    r"\label{tab:rates}"
    "\n"
    r"\begin{tabular}{ccrrr}"
    "\n"
    r"\toprule"
    "\n"
    r"Class $k$ & $c_k$ & AIC (exp) & AIC (power) & Preferred \\"
    "\n"
    r"\midrule"
    "\n"
    r"0 & 0.750 & 45.5 & 79.1 & Exponential \\"
    "\n"
    r"1 & 0.985 & 55.8 & 104.8 & Exponential \\"
    "\n"
    r"2 & 0.401 & 6.2 & 10.1 & Exponential \\"
    "\n"
    r"3 & 0.401 & 6.2 & 10.1 & Exponential \\"
    "\n"
    r"4 & 0.985 & 55.8 & 104.8 & Exponential \\"
    "\n"
    r"\bottomrule"
    "\n"
    r"\end{tabular}"
    "\n"
    r"\end{table}"
    "\n\n"
    r"This exponential approach to zero is expected from generic hyperbolic dynamics"
    "\n"
    r"(closed geodesics become dense), but the class-dependent rates encode"
    "\n"
    r"the torsion structure."
)

tex = tex.replace(old_decay_eq, new_decay_eq, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("CMP format applied, conjecture and rate table added.")
