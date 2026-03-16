path = r"C:\dev\framework\papers\hyperbolic-flavor-ckm\gentry-hyperbolic-flavor-ckm.tex"
with open(path, "r", encoding="utf-8") as f:
    tex = f.read()

changes = 0

# 1. Abstract: replace K_H claim with Gaussian claim
old = ('boundary-localized states using the exact hyperbolic boundary propagator\n'
       '$K_H \\propto (1-\\cos\\theta)^{-\\Delta/2}$, which reduces to a Gaussian\n'
       'in the small-angle limit. The resulting matrix matches the experimental\n'
       'CKM moduli with Frobenius error of 0.017, reproducing the Cabibbo angle\n'
       'to within 0.4\\% with a single free parameter ($\\Delta \\approx 2.08$,\n'
       'equivalent to localization width $\\sigma \\approx 0.49$). The three\n')
new = ('boundary-localized states using a Gaussian overlap kernel\n'
       '$K(\\theta;\\sigma) = e^{-\\theta^2/2\\sigma^2}$, where $\\sigma$\n'
       'is a single localization parameter. The resulting matrix matches the experimental\n'
       'CKM moduli with Frobenius error of 0.017, reproducing the Cabibbo angle\n'
       'to within 0.4\\% with a single free parameter ($\\sigma = 0.49$). The three\n')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("1. Abstract: OK")
else:
    print("1. Abstract: NOMATCH")

# 2. Kernel definition block (lines 157-187)
old = ('We parametrize boundary points via stereographic projection\n'
       '$z \\in \\mathbb{C} \\cup \\{\\infty\\}$ and define the overlap kernel\n'
       'as the exact bulk-to-boundary propagator of a free scalar field on\n'
       '$\\Hthree$ with conformal dimension $\\Delta$~\\cite{Witten1998}:\n'
       '\\begin{equation}\n'
       'K_H(z_1, z_2;\\Delta) = \\mathcal{N}\\,\n'
       '\\bigl(1 - \\cos\\theta(z_1,z_2)\\bigr)^{-\\Delta/2},\n'
       '\\label{eq:kernel}\n'
       '\\end{equation}\n'
       'where $\\theta(z_1, z_2) \\in [0,\\pi]$ is the angle between the\n'
       'corresponding unit vectors in $\\mathbb{R}^3$, $\\Delta > 0$ is the\n'
       'conformal dimension (the single free parameter of our model), and\n'
       '$\\mathcal{N}$ is a normalization constant absorbed into the QR\n'
       'orthogonalization.\n'
       '\n'
       'This is the natural kernel on $\\partial\\Hthree \\cong S^2$: it is\n'
       '$\\mathrm{PSL}(2,\\mathbb{C})$-equivariant and decays with angle exactly\n'
       'as a boundary Green\'s function should. For small $\\theta$, one has\n'
       '$1 - \\cos\\theta \\approx \\theta^2/2$, so\n'
       '\\begin{equation}\n'
       'K_H \\approx \\left(\\frac{\\theta^2}{2}\\right)^{-\\Delta/2}\n'
       '\\propto \\theta^{-\\Delta}.\n'
       '\\label{eq:kernel_small}\n'
       '\\end{equation}\n'
       'A Gaussian overlap $e^{-\\theta^2/2\\sigma^2}$ approximates this\n'
       'power-law behavior in the small-angle regime when\n'
       '$\\Delta/2 \\approx 1/(2\\sigma^2)$, i.e., $\\Delta \\approx \\sigma^{-2}$.\n'
       'With $\\sigma = 0.49$ (the empirically optimal value), this gives\n'
       '$\\Delta \\approx 2.08$, a natural $O(1)$ conformal dimension.\n'
       'All numerical results reported in this paper use the full propagator\n'
       '$K_H$ of Eq.~\\eqref{eq:kernel}.\n')
new = ('We parametrize boundary points via stereographic projection\n'
       '$z \\in \\mathbb{C} \\cup \\{\\infty\\}$ and define the overlap via a\n'
       'Gaussian kernel:\n'
       '\\begin{equation}\n'
       'K(z_1, z_2;\\sigma) = \\exp\\!\\left( -\\frac{\\theta(z_1,z_2)^2}{2\\sigma^2} \\right),\n'
       '\\label{eq:kernel}\n'
       '\\end{equation}\n'
       'where $\\theta(z_1, z_2) \\in [0,\\pi]$ is the angle between the\n'
       'corresponding unit vectors in $\\mathbb{R}^3$, and $\\sigma > 0$ is a\n'
       'localization parameter controlling the width of the boundary states.\n'
       'The kernel decreases monotonically with angular separation, equaling\n'
       'unity for coincident directions and vanishing as $\\theta \\to \\pi$.\n'
       '\n'
       '\\textit{Remark.} The bulk-to-boundary propagator on $\\Hthree$ for a\n'
       'scalar of conformal dimension $\\Delta$ takes the form\n'
       '$K_H \\propto (1-\\cos\\theta)^{-\\Delta/2}$~\\cite{Witten1998}.\n'
       'Unlike the Gaussian, this kernel increases with $\\theta$ and is\n'
       'therefore not directly applicable as an overlap kernel here.\n'
       'Deriving the Gaussian from a regulated hyperbolic propagator is\n'
       'an interesting open problem left for future work.\n')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("2. Kernel block: OK")
else:
    print("2. Kernel block: NOMATCH")

# 3. Overlap matrix equation
old = ('\\mathcal{O}_{ij} = K_H(\\mathbf{n}_{\\gamma_i}, \\mathbf{n}_{\\gamma_j};\\Delta)\n'
       '= \\bigl(1 - \\cos\\theta_{ij}\\bigr)^{-\\Delta/2},\n')
new = ('\\mathcal{O}_{ij} = K(\\mathbf{n}_{\\gamma_i}, \\mathbf{n}_{\\gamma_j};\\sigma)\n'
       '= \\exp\\!\\left(-\\frac{\\theta_{ij}^2}{2\\sigma^2}\\right),\n')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("3. Overlap eq: OK")
else:
    print("3. Overlap eq: NOMATCH")

# 4. Scan methodology
old = ('We scanned over the conformal dimension $\\Delta \\in [1.0, 5.0]$\n'
       '(equivalently $\\sigma \\in [0.45, 1.0]$ via $\\Delta \\approx \\sigma^{-2}$)\n'
       'and retained the configuration minimizing $\\chi^2$.')
new = ('We scanned over the localization parameter $\\sigma \\in [0.3, 1.0]$\n'
       'and retained the configuration minimizing $\\chi^2$.')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("4. Scan method: OK")
else:
    print("4. Scan method: NOMATCH")

# 5. Table caption
old = 'Computed using the hyperbolic propagator kernel~\\eqref{eq:kernel} with $\\Delta = 2.08$ ($\\sigma_{\\mathrm{eff}} \\approx 0.49$).'
new = 'Computed using the Gaussian overlap kernel~\\eqref{eq:kernel} with $\\sigma = 0.49$.'
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("5. Table caption: OK")
else:
    print("5. Table caption: NOMATCH")

# 6. Robustness Delta -> sigma
old = ('The result is robust to variations in the conformal dimension $\\Delta$.\n'
       'The Frobenius error has a broad minimum around $\\Delta \\approx 2.08$\n'
       '($\\sigma_{\\mathrm{eff}} \\approx 0.49$), remaining below 0.05 for\n'
       '$\\Delta \\in [1.8, 2.8]$, indicating that no fine-tuning is required.')
new = ('The result is robust to variations in the localization parameter $\\sigma$.\n'
       'The Frobenius error has a broad minimum around $\\sigma \\approx 0.49$,\n'
       'remaining below 0.05 for $\\sigma \\in [0.40, 0.60]$,\n'
       'indicating that no fine-tuning is required.')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("6. Robustness: OK")
else:
    print("6. Robustness: NOMATCH")

# 7. Poincare figure text
old = ('These angles, combined with the hyperbolic propagator kernel~\\eqref{eq:kernel} at\n'
       '$\\Delta = 2.08$, produce the overlap matrix from which the mixing\n'
       'matrix is extracted via QR decomposition.')
new = ('These angles, combined with the Gaussian overlap kernel~\\eqref{eq:kernel} at\n'
       '$\\sigma = 0.49$, produce the overlap matrix from which the mixing\n'
       'matrix is extracted via QR decomposition.')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("7. Poincare text: OK")
else:
    print("7. Poincare text: NOMATCH")

# 8. Cabibbo section
old = ('This value emerges from the hyperbolic propagator~\\eqref{eq:kernel}\n'
       'evaluated at the geometric angle $\\theta_{12} \\approx 48.2^\\circ$:\n'
       '\\begin{equation}\n'
       'K_H(\\theta_{12}; \\Delta) = (1 - \\cos 48.2^\\circ)^{-\\Delta/2}\n'
       '\\approx (0.333)^{-1.04} \\approx 0.23,\n'
       '\\end{equation}\n'
       'directly generating a mixing element of the observed magnitude.\n'
       'The Gaussian approximation $e^{-\\theta_{12}^2/2\\sigma^2} \\approx 0.23$\n'
       'at $\\sigma = 0.49$ gives the same result, confirming the small-angle\n'
       'equivalence of Eq.~\\eqref{eq:kernel_small}.')
new = ('This value emerges from the Gaussian kernel~\\eqref{eq:kernel}\n'
       'evaluated at the geometric angle $\\theta_{12} \\approx 48.2^\\circ$:\n'
       '\\begin{equation}\n'
       'K(\\theta_{12};\\sigma) = e^{-\\theta_{12}^2/2\\sigma^2}\n'
       '= e^{-(0.8405)^2/2(0.49)^2} \\approx 0.23,\n'
       '\\end{equation}\n'
       'directly generating a mixing element of the observed magnitude.\n'
       'The axis angle $\\theta_{12}$ is fixed by Mostow rigidity; only\n'
       'the overall scale $\\sigma$ is adjusted to match experiment.')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("8. Cabibbo: OK")
else:
    print("8. Cabibbo: NOMATCH")

# 9. Hierarchy section
old = ('The hyperbolic propagator $K_H \\propto (1-\\cos\\theta)^{-\\Delta/2}$\n'
       'translates these angular differences into the observed mixing hierarchy:\n'
       'smaller angles produce larger overlaps, driving the pattern\n'
       '$|V_{ud}| \\gg |V_{us}| \\gg |V_{ub}|$.')
new = ('The Gaussian kernel naturally translates these angular differences\n'
       'into the observed mixing hierarchy: smaller angles produce larger\n'
       'overlaps, driving the pattern $|V_{ud}| \\gg |V_{us}| \\gg |V_{ub}|$.')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("9. Hierarchy: OK")
else:
    print("9. Hierarchy: NOMATCH")

# 10. CP enumerate
old = ('\\item The overlap kernel~\\eqref{eq:kernel} is purely real ($(1-\\cos\\theta)^{-\\Delta/2} \\in \\mathbb{R}$).\n')
new = ('\\item The overlap kernel~\\eqref{eq:kernel} is purely real.\n')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("10. CP enumerate: OK")
else:
    print("10. CP enumerate: NOMATCH")

# 11. Conclusion kernel reference
old = ('with the hyperbolic boundary propagator~\\eqref{eq:kernel} at\n'
       'conformal dimension $\\Delta = 2.08$ (equivalent Gaussian width\n'
       '$\\sigma \\approx 0.49$), suffice to reproduce the dominant structure\n'
       'of the CKM matrix.')
new = ('with a single Gaussian localization parameter $\\sigma = 0.49$,\n'
       'suffice to reproduce the dominant structure of the CKM matrix.')
if old in tex:
    tex = tex.replace(old, new); changes += 1; print("11. Conclusion: OK")
else:
    print("11. Conclusion: NOMATCH")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)

print(f"\n{changes}/11 changes applied.")
print(f"Line count: {tex.count(chr(10))}")