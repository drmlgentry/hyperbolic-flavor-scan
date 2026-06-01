path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()
changes = []

# Fix 1: Abstract MZ/MW wording
if r'appears at $0.09\%$' in tex:
    tex = tex.replace(r'appears at $0.09\%$', r'is reproduced to $0.09\%$')
    changes.append('Fix 1: abstract MZ/MW wording')

# Fix 2: Introduction "companion series" -> "previous work"
if r'In a companion series~\cite{GentryCKM,GentryPMNS,GentryCP}' in tex:
    tex = tex.replace(
        r'In a companion series~\cite{GentryCKM,GentryPMNS,GentryCP}',
        r'In previous work~\cite{GentryCKM,GentryPMNS,GentryCP}'
    )
    changes.append('Fix 2: intro wording')
elif r'In related work submitted to Physical Review' in tex:
    tex = tex.replace(
        r'In related work submitted to Physical Review',
        r'In previous work submitted to Physical Review'
    )
    changes.append('Fix 2: intro wording (alternate)')

# Fix 3: Table 1 MZ/MW error -- add absolute value
tex = tex.replace(
    r'$M_Z/M_W=1.13451$ & \mfld{m006} & \word{aaabb/baaab} & $1.1355$ & r & $0.09\%$ \\',
    r'$M_Z/M_W=1.13451$ & \mfld{m006} & \word{aaabb/baaab} & $1.1355$ & r & $0.0010$ (0.09\%) \\'
)
changes.append('Fix 3: Table 1 MZ/MW absolute error')

# Fix 4: Add forward reference to appendix in Selberg section
old_selberg = r"""\textbf{Convergence.}
The truncation introduces systematic errors that we have not rigorously
bounded."""
new_selberg = r"""\textbf{Convergence.}
Full details of the numerical method, truncation criteria, and remainder
bounds are given in Appendix~\ref{app:selberg}.
The truncation introduces systematic errors that we have not rigorously
bounded."""
if old_selberg in tex:
    tex = tex.replace(old_selberg, new_selberg)
    changes.append('Fix 4: forward ref to appendix added in Selberg section')

# Fix 5: Add uncertainty qualifier to lambda1 results
old_lam = r"""\textbf{Results.}
\begin{align}
  \lambda_1(\mfld{m003}) &\approx 2.48 \quad (r_1\approx 1.49),
  \label{eq:lam1_m003}\\
  \lambda_1(\mfld{m006}) &\approx 2.82 \quad (r_1\approx 1.60).
  \label{eq:lam1_m006}
\end{align}"""
new_lam = r"""\textbf{Results.}
These are preliminary estimates based on a truncated product;
as detailed in Appendix~\ref{app:selberg}, the values are stable
to approximately $\pm 0.1$ under variation of the truncation parameters.
\begin{align}
  \lambda_1(\mfld{m003}) &\approx 2.48\pm 0.1 \quad (r_1\approx 1.49),
  \label{eq:lam1_m003}\\
  \lambda_1(\mfld{m006}) &\approx 2.82\pm 0.1 \quad (r_1\approx 1.60).
  \label{eq:lam1_m006}
\end{align}"""
if old_lam in tex:
    tex = tex.replace(old_lam, new_lam)
    changes.append('Fix 5: lambda1 uncertainty qualifier added')

for c in changes: print(c)
with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
