path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

old = r'''\subsection{Numerical results}

Table~\ref{tab:complex_results} summarizes the best results
from the complex Borel construction scan on $m003$.

\begin{table}[h]
\centering
\caption{Best results from complex Borel construction on $m003$.
$\alpha$ is the imaginary weight~\eqref{eq:complex_axis},
$\mathcal{F}$ is permutation-minimized Frobenius fitness,
$J$ is the Jarlskog invariant.
PDG targets: $\mathcal{F} = 0.01897$, $J = -0.00966$.}
\label{tab:complex_results}
\begin{tabular}{llccc}
\toprule
Words & $\alpha$ & scale & $\mathcal{F}$ & $J$ \\
\midrule
aa/ab/aB (real) & $0$ & $1.0$ & $0.01897$ & $0.000$ \\
Ab/aB/aa        & $0.897$ & $1.729$ & $0.0489$ & $-0.017$ \\
Ab/aa/aB        & $1.198$ & $1.916$ & $0.0931$ & $+0.020$ \\
\bottomrule
\end{tabular}
\end{table}

The best result (words \texttt{Ab/aB/aa}, $\alpha=0.897$)
achieves $J = -0.017$, which has the correct sign and is within
a factor of $1.75$ of the target $J_{\mathrm{PDG}} = -0.00966$.
The fitness degradation from $0.019$ to $0.049$ reflects a
genuine tension: the real construction optimizes $|U_{\mathrm{PMNS}}|$
at $\alpha=0$ while the complex construction optimizes $J \neq 0$
at $\alpha \neq 0$.'''

new = r'''\subsection{Numerical results}

Table~\ref{tab:complex_results} summarizes the best results
from the complex Borel construction census scan on $m003$,
using mixed word triples (two length-2 words and one length-3 word).
The scan covered 79,872 triples across six values of $\alpha$.

\begin{table}[h]
\centering
\caption{Best results from complex Borel construction on $m003$
with length-3 words. $\alpha$ is the imaginary weight~\eqref{eq:complex_axis}.
PDG targets: $\mathcal{F} = 0.01897$, $J_{\mathrm{PDG}} = -0.00966$.}
\label{tab:complex_results}
\begin{tabular}{llccc}
\toprule
Words & $\alpha$ & scale & $\mathcal{F}$ & $J$ \\
\midrule
aa/ab/aB (real)   & $0$     & $1.000$ & $0.01897$ & $0.000$ \\
Ba/baB/BBa        & $0.900$ & $2.095$ & $0.01990$ & $-0.01256$ \\
Ba/baB/BBa        & $1.100$ & $2.095$ & $0.03855$ & $-0.00977$ \\
Ba/baB/BBa        & $1.125$ & $2.132$ & $0.03689$ & $-0.00966$ \\
\bottomrule
\end{tabular}
\end{table}

The optimal word triple \texttt{Ba/baB/BBa} (equivalently
\texttt{Ab/bAB/Abb}) achieves a Pareto-optimal trade-off
between fitness and $J$:
\begin{itemize}
\item At $\alpha = 0.900$, scale $= 2.095$: fitness $= 0.0199$,
matching the real construction, with $J = -0.01256$
(within a factor of $1.3$ of $J_{\mathrm{PDG}}$).
\item At $\alpha = 1.125$, scale $= 2.132$: $J = -0.00966$,
matching $J_{\mathrm{PDG}}$ to four significant figures,
with fitness $= 0.037$.
\end{itemize}
The output matrix at the fitness-optimal point ($\alpha = 0.900$) is
\begin{equation}
|U_{\mathrm{geom}}^{(\mathrm{c})}| = \begin{pmatrix}
0.8201 & 0.5529 & 0.1473 \\
0.3610 & 0.3233 & 0.8747 \\
0.4439 & 0.7680 & 0.4617
\end{pmatrix},
\end{equation}
compared to the PDG values~\eqref{eq:U_pmns_pdg}.
This demonstrates that a single geometric construction with
$\alpha = 0.9$ and length-3 words simultaneously reproduces
both $|U_{\mathrm{PMNS}}|$ to within $\mathcal{F} = 0.020$ and
generates $J$ with the correct sign and order of magnitude.
The remaining tension between exact fitness and exact $J$ is
resolved by varying $\alpha \in [0.9, 1.125]$, tracing a
Pareto frontier in the $(\mathcal{F}, J)$ plane.'''

if old in tex:
    tex = tex.replace(old, new)
    print('Section 6.3 updated')
else:
    print('Pattern not found')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
