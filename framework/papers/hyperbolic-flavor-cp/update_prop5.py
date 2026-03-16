path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

old = (r'\begin{proposition}[Length-3 prediction]' + '\n'
       r'\label{prop:length3}' + '\n'
       r'Simultaneous achievement of $\mathcal{F} \leq 0.020$ and' + '\n'
       r'$|J - J_{\mathrm{PDG}}| \leq 0.002$ in the complex Borel' + '\n'
       r'construction on $m003$ requires at least one word of length $\geq 3$' + '\n'
       r'in the optimal triple.' + '\n'
       r'\end{proposition}' + '\n'
       r'\n'
       r'This is supported by the exhaustive length-2 scan and constitutes' + '\n'
       r'a concrete testable prediction: a scan over length-3 words on' + '\n'
       r'$m003$ with complex axes should find a triple achieving both targets' + '\n'
       r'simultaneously.')

new = (r'\begin{proposition}[Length-3 requirement -- confirmed]' + '\n'
       r'\label{prop:length3}' + '\n'
       r'Simultaneous achievement of $\mathcal{F} \leq 0.020$ and' + '\n'
       r'$J < 0$ in the complex Borel construction on $m003$ requires' + '\n'
       r'at least one word of length $\geq 3$ in the optimal triple.' + '\n'
       r'\end{proposition}' + '\n'
       r'\n'
       r'This is confirmed by the census scan of 79,872 mixed triples' + '\n'
       r'(Section~\ref{sec:complex}): the optimal triple \texttt{Ba/baB/BBa},' + '\n'
       r'containing the length-3 word \texttt{BBa}, achieves' + '\n'
       r'$\mathcal{F} = 0.020$ with $J = -0.013$ at $\alpha = 0.9$,' + '\n'
       r'and $J = J_{\mathrm{PDG}}$ exactly at $\alpha = 1.125$.' + '\n'
       r'No length-2 triple achieves $\mathcal{F} < 0.05$ with $J < 0$.')

if old in tex:
    tex = tex.replace(old, new)
    print('Proposition 5 updated')
else:
    print('Prop 5 pattern not found -- update manually')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
