path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# (h) abstract: "the" row-phase
tex = tex.replace(
    'identify the obstruction as a row-phase invariance theorem,',
    'identify the obstruction as the row-phase invariance theorem,')
print('(h) done' if 'the row-phase invariance theorem' in tex else '(h) NOT FOUND')

# (a) alpha>=0 and n_Im clarification
tex = tex.replace(
    r'where $\alpha \in [0,1]$ interpolates between the real ($\alpha=0$)' + '\n' +
    r'and fully complex ($\alpha=1$) constructions, and' + '\n' +
    r'$\hat{n}^{(\mathrm{Im})}$ is the normalized imaginary part.',
    r'where $\alpha \geq 0$ is the imaginary weight (values $\alpha > 1$ are' + '\n' +
    r'permitted; the denominator remains well-defined for all $\alpha \geq 0$),' + '\n' +
    r'$\hat{n}^{(\mathrm{Re})}$ is the unit vector from the real Pauli components' + '\n' +
    r'of $\log\rho(w)$ (as in the real construction), and' + '\n' +
    r'$\hat{n}^{(\mathrm{Im})}$ is the analogous unit vector from the imaginary' + '\n' +
    r'Pauli components $(\mathrm{Im}[L_{01}+L_{10}],$' + '\n' +
    r'$\mathrm{Im}[i(L_{10}-L_{01})], \mathrm{Im}[L_{00}-L_{11}])$, normalized.')
print('(a) done' if r'\alpha \geq 0' in tex else '(a) NOT FOUND')

# (e) isospectrality elaboration
tex = tex.replace(
    'reflecting the arithmetic structure of its fundamental group.',
    'reflecting the arithmetic structure of its fundamental group.\n' +
    r'Geometrically, the geodesics for $a$ and $aB = a \cdot b^{-1}$ have' + '\n' +
    r'equal complex length; this may reflect a symmetry of $m003$ that' + '\n' +
    r'exchanges $b$ and its inverse $B$. The precise group-theoretic' + '\n' +
    r'reason is an open question.')
print('(e) done' if r'exchanges $b$ and its inverse' in tex else '(e) NOT FOUND')

# (f) Pareto figure caption: alpha>1 note
tex = tex.replace(
    'Color encodes $\\alpha$.',
    'Color encodes $\\alpha$. Values $\\alpha > 1$ are geometrically\n' +
    r'permitted since the norm in~\eqref{eq:complex_axis} is well-defined' + '\n' +
    r'for all $\alpha \geq 0$.')
print('(f) done' if r'Values $\\alpha > 1$' in tex else '(f) NOT FOUND')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('All done.')
