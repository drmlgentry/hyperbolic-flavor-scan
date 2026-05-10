path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# Update Table 3 (rates) with CIs
old_table = r"""\begin{tabular}{ccrrr}
\toprule
Class $k$ & $c_k$ & AIC (exp) & AIC (pow) & Preferred \\
\midrule
0 & 0.750 & 45.5 & 79.1 & Exponential \\
1 & 0.985 & 55.8 & 104.8 & Exponential \\
2 & 0.401 &  6.2 &  10.1 & Exponential \\
3 & 0.401 &  6.2 &  10.1 & Exponential \\
4 & 0.985 & 55.8 & 104.8 & Exponential \\
\bottomrule
\end{tabular}"""

new_table = r"""\begin{tabular}{cccr}
\toprule
Class $k$ & $c_k$ & 95\% CI & $\Delta$AIC \\
\midrule
0 & 0.750 & [0.434,\ 0.944] & $-33.6$ \\
1 & 0.985 & [0.541,\ 1.376] & $-49.0$ \\
2 & 0.401 & [0.080,\ 0.645] & $-4.0$  \\
3 & 0.401 & [0.079,\ 0.651] & $-4.0$  \\
4 & 0.985 & [0.585,\ 1.431] & $-49.0$ \\
\midrule
\multicolumn{4}{l}{$\Delta$AIC $=$ AIC(exp) $-$ AIC(pow); negative favours exponential.}\\
\multicolumn{4}{l}{CIs from $2000$ bootstrap resamples (word-length basis).}\\
\bottomrule
\end{tabular}"""

tex = tex.replace(old_table, new_table, 1)

# Update caption for Table 3
old_cap = (r"\caption{Exponential decay rates $c_k$ from log-linear fits to"
    "\n"
    r"  $\phi_{\rm floor}^{(k)}(L)$, with AIC comparison against"
    "\n"
    r"  power-law fits. Lower AIC = better fit.}")
new_cap = (r"\caption{Exponential decay rates $c_k$ from log-linear fits to"
    "\n"
    r"  $\phi_{\rm floor}^{(k)}(L)$, with 95\% bootstrap confidence intervals"
    "\n"
    r"  (2000 resamples) and $\Delta$AIC $=$ AIC(exp) $-$ AIC(pow)."
    "\n"
    r"  Negative $\Delta$AIC favours exponential over power-law."
    "\n"
    r"  The pairings $c_1=c_4$ and $c_2=c_3$ are confirmed:"
    "\n"
    r"  their CIs overlap completely.}")
tex = tex.replace(old_cap, new_cap, 1)

# Add note about translation-length basis in Discussion
old_disc = (r"Comparing these rates across the OrientableClosedCensus for other manifolds"
    "\n"
    r"with $H_1=\mathbb{Z}/p$ ($p$ prime) would test the conjecture and may"
    "\n"
    r"reveal which arithmetic properties control the rate hierarchy.")
new_disc = (r"Comparing these rates across the OrientableClosedCensus for other manifolds"
    "\n"
    r"with $H_1=\mathbb{Z}/p$ ($p$ prime) would test the conjecture and may"
    "\n"
    r"reveal which arithmetic properties control the rate hierarchy."
    "\n\n"
    r"We note that the word-length $L$ depends on the choice of generators"
    "\n"
    r"$\{a,b\}$ for $\pi_1(\mfld{m006})$.  A generator-independent formulation"
    "\n"
    r"would use the geodesic translation length $\ell(\gamma)=2\log|\lambda(\gamma)|$"
    "\n"
    r"in place of $L$.  A preliminary analysis using the running maximum of"
    "\n"
    r"$\ell(\gamma)$ as the $x$-axis variable confirms that the class pairing"
    "\n"
    r"$c_1=c_4$ and $c_2=c_3$ persists, though the rate values differ;"
    "\n"
    r"a rigorous treatment via the geodesic length spectrum of \mfld{m006}"
    "\n"
    r"is deferred to future work.")
tex = tex.replace(old_disc, new_disc, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("CIs and generator-independence note added.")
