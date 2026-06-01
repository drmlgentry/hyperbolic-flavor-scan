path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: lines = f.readlines()

replacements = {
    43:  r'$\delta_{\mathrm{CP}} = \phi_{aa} - \phi_{ab} + \phi_{aB} = 203.5^\circ$,' + '\n',
    44:  r'within $3.3\%$ of the PDG 2024 value of $197^\circ$ with no free parameters.' + '\n',
    93:  r'The combination $\phi_{aa} - \phi_{ab} + \phi_{aB} = 203.5^\circ$' + '\n',
    160: r'\texttt{aa} & $-0.8894$ & $-168.6^\circ$ & $2.4338$ \\' + '\n',
    161: r'\texttt{ab} & $-0.4992$ & $+83.7^\circ$ & $1.6473$ \\' + '\n',
    162: r'\texttt{aB} & $-0.4447$ & $+95.7^\circ$ & $1.5601$ \\' + '\n',
    194: r'\delta_{\mathrm{geom}} &= (-168.6^\circ) - (83.7^\circ) + (95.7^\circ) \nonumber \\' + '\n',
    195: r'&= -156.6^\circ \equiv 203.5^\circ \pmod{360^\circ}.' + '\n',
    198: r'The PDG 2024 central value is $\delta_{\mathrm{CP}} = 197^\circ$,' + '\n',
    199: r'giving a discrepancy of $6.5^\circ$ or $3.3\%$.' + '\n',
    213: r'$\phi_{aa}/2 = \phi_a = -84.3^\circ$ is the fundamental twist of' + '\n',
    421: r'\texttt{aaB} & $+0.8859$ & $92.49^\circ$ & $2.425$ \\' + '\n',
    422: r'\texttt{AbA} & $+0.8859$ & $92.49^\circ$ & $2.425$ \\' + '\n',
    423: r'\texttt{AAb} & $+0.8859$ & $92.49^\circ$ & $2.425$ \\' + '\n',
    432: r'$t = 0.8859$ and $\phi = 92.49^\circ \approx \pi/2$.' + '\n',
    458: r'The common twist angle $\phi \approx 92.49^\circ \approx \pi/2$' + '\n',
    461: r'$7\theta_{12}^{\mathrm{CKM}} = 7 \times 13.04^\circ = 91.28^\circ$,' + '\n',
    462: r'within $1.2^\circ$ of $\phi$, suggesting a possible connection' + '\n',
    486: r'$\delta_{\mathrm{CP}} = 203.5^\circ$, within $3.3\%$ of the PDG value' + '\n',
    487: r'of $197^\circ$, with no free parameters.' + '\n',
    583: r'All three words are isospectral with $\phi \approx 92.5^\circ$.}' + '\n',
}

for lineno, newline in replacements.items():
    idx = lineno - 1
    lines[idx] = newline
    print('Fixed line ' + str(lineno))

with open(path, 'w', encoding='utf-8') as f: f.writelines(lines)
print('All done.')
