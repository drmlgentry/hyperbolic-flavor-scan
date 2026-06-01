path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

fixes = [

# (a) Clarify complex axis definition -- add sentence after eq:complex_axis
(
r'where $\alpha \in [0,1]$ interpolates between the real ($\alpha=0$)'
r'\nand fully complex ($\alpha=1$) constructions, and'
r'\n$\hat{n}^{(\mathrm{Im})}$ is the normalized imaginary part.',

r'where $\alpha \geq 0$ is the imaginary weight (values $\alpha > 1$ are'
r'\npermitted; the denominator remains well-defined for all $\alpha$),'
r'\n$\hat{n}^{(\mathrm{Re})}$ is the unit vector from the real Pauli'
r'\ncomponents of $\log\rho(w)$ (as in the real construction), and'
r'\n$\hat{n}^{(\mathrm{Im})}$ is the unit vector from the imaginary Pauli'
r'\ncomponents: $\mathbf{n}^{(\mathrm{Im})} = (\mathrm{Im}[L_{01}+L_{10}],'
r'\n\mathrm{Im}[i(L_{10}-L_{01})], \mathrm{Im}[L_{00}-L_{11}])$,'
r'\nnormalized to unit length.'
),

# (b) Add geometric explanation for CP combination
(
r'The combination $\phi_{aa} - \phi_{ab} + \phi_{aB}$ has a natural'
r'\ngeometric interpretation. The word $aa = a^2$ satisfies'
r'\n$\phi(a^2) = 2\phi(a)$ by additivity of twist angles'
r'\n(Proposition~\ref{prop:length_theorem}), so'
r'\n$\phi_{aa}/2 = \phi_a = -84.3^\circ$ is the fundamental twist of'
r'\nthe generator $a$. The combination becomes:'
r'\n\begin{equation}'
r'\n\delta_{\mathrm{geom}} = 2\phi_a - \phi_{ab} + \phi_{aB},'
r'\n\end{equation}'
r'\nwhich measures the signed deficit of twist angles around the'
r'\ngenerator triangle $(a, b, B)$ on $m003$.'
r'\nThe positive sign for $\phi_{aa}$ (the squared generator) and'
r'\nmixed signs for the cross-words reflect the orientation of the'
r'\ngeodesic triangle: traversing $a \to ab \to aB \to a$ accumulates'
r'\na net twist equal to $\delta_{\mathrm{geom}}$.',

r'The combination $\phi_{aa} - \phi_{ab} + \phi_{aB}$ has a natural'
r'\ngeometric interpretation. The word $aa = a^2$ satisfies'
r'\n$\phi(a^2) = 2\phi(a)$ by additivity of twist angles'
r'\n(Proposition~\ref{prop:length_theorem}), so'
r'\n$\phi_{aa}/2 = \phi_a = -84.3^\circ$ is the fundamental twist of'
r'\nthe generator $a$. The signed combination'
r'\n$\phi_{aa} - \phi_{ab} + \phi_{aB} = 2\phi_a - \phi_{ab} + \phi_{aB}$'
r'\nmeasures the net twist accumulated by traversing the geodesic triangle'
r'\n$(a, ab, aB)$ on $m003$: the positive weight on $\phi_{aa}$'
r'\nreflects the squared generator, while the alternating signs on'
r'\n$\phi_{ab}$ and $\phi_{aB}$ reflect the two cross-words.'
r'\nBy the additivity of twist angles along concatenated geodesics,'
r'\nthis combination is the unique rephasing-invariant combination'
r'\nconsistent with the orientation of the triangle.'
),

# (c) Clarify Table 1 caption -- eigenvalue magnitude
(
r'$e^{|t|}$ is the eigenvalue magnitude $|\lambda| = e^{|t|}$.',
r'$|\lambda| = e^{|t|}$ is the eigenvalue modulus.'
),

# (e) Add note on isospectrality of a and aB
(
r'This is a non-trivial constraint on the length spectrum of $m003$,'
r'\nreflecting the arithmetic structure of its fundamental group.',

r'This is a non-trivial constraint on the length spectrum of $m003$.'
r'\nGeometrically, it means that the geodesics corresponding to $a$'
r'\nand $aB = a \cdot b^{-1}$ have the same complex length.'
r'\nThis may reflect a symmetry of $m003$ that exchanges the roles'
r'\nof $b$ and its inverse $B$, or equivalently that the holonomy'
r'\nof $b^{-1}$ contributes zero net translation length when composed'
r'\nwith $a$. The precise group-theoretic reason is an open question.'
),

# (f) Add alpha>1 note to Pareto figure caption
(
r'Color encodes $\\alpha$.',
r'Color encodes $\\alpha$. Values $\\alpha > 1$ are geometrically'
r'\npermitted since the complex axis norm~\\eqref{eq:complex_axis}'
r'\nremains well-defined for all $\\alpha \\geq 0$.'
),

# (h) Abstract: "the" row-phase invariance theorem
(
r'identify the obstruction as a row-phase invariance theorem,',
r'identify the obstruction as the row-phase invariance theorem,'
),

]

count = 0
for old, new in fixes:
    old_clean = old.replace(r'\n', '\n')
    new_clean = new.replace(r'\n', '\n')
    if old_clean in tex:
        tex = tex.replace(old_clean, new_clean)
        count += 1
        print('Fixed: ' + old[:60].replace('\n',' '))
    else:
        print('NOT FOUND: ' + old[:60].replace('\n',' '))

# (c) Table 1 caption fix -- find and update
import re
tex = re.sub(
    r'\$e\^\{\\?|t\\?\}\$\s+is the eigenvalue magnitude',
    r'$|\\lambda| = e^{|t|}$ is the eigenvalue modulus',
    tex
)

print()
print(str(count) + ' fixes applied')
with open(path, 'w', encoding='utf-8') as f: f.write(tex)
