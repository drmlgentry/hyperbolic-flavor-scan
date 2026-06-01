path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Fix the PMNS A-factor table rows
# Current broken: \texttt{aa} & $-0.8894$ & $-$168.6\ensuremath{^\circ}$ & $2.4338$ \\
# Correct:        \texttt{aa} & $-0.8894$ & $-168.6^\circ$ & $2.4338$ \\
fixes = [
    # PMNS table
    (r'\texttt{aa} & $-0.8894$ & $-$168.6\ensuremath{^\circ}$ & $2.4338$ \\',
     r'\texttt{aa} & $-0.8894$ & $-168.6^\circ$ & $2.4338$ \\'),
    (r'\texttt{ab} & $-0.4992$ & $+$83.7\ensuremath{^\circ}$  & $1.6473$ \\',
     r'\texttt{ab} & $-0.4992$ & $+83.7^\circ$ & $1.6473$ \\'),
    (r'\texttt{aB} & $-0.4447$ & $+$95.7\ensuremath{^\circ}$  & $1.5601$ \\',
     r'\texttt{aB} & $-0.4447$ & $+95.7^\circ$ & $1.5601$ \\'),
    # CKM table
    (r'\texttt{aaB} & $+0.8859$ & $$92.49\ensuremath{^\circ}$ & $2.425$ \\',
     r'\texttt{aaB} & $+0.8859$ & $92.49^\circ$ & $2.425$ \\'),
    (r'\texttt{AbA} & $+0.8859$ & $$92.49\ensuremath{^\circ}$ & $2.425$ \\',
     r'\texttt{AbA} & $+0.8859$ & $92.49^\circ$ & $2.425$ \\'),
    (r'\texttt{AAb} & $+0.8859$ & $$92.49\ensuremath{^\circ}$ & $2.425$ \\',
     r'\texttt{AAb} & $+0.8859$ & $92.49^\circ$ & $2.425$ \\'),
    # Appendix table caption
    (r'isospectral with $\phi \approx $92.5\ensuremath{^\circ}$.}',
     r'isospectral with $\phi \approx 92.5^\circ$.}'),
    # Inline text occurrences
    (r'$\phi = $92.49\ensuremath{^\circ}$ \approx \pi/2$',
     r'$\phi = 92.49^\circ \approx \pi/2$'),
    (r'$\phi \approx $92.49\ensuremath{^\circ}$ \approx \pi/2$',
     r'$\phi \approx 92.49^\circ \approx \pi/2$'),
    (r'$7 \times $13.04\ensuremath{^\circ}$ = $91.28\ensuremath{^\circ}$$',
     r'$7 \times 13.04^\circ = 91.28^\circ$'),
    (r'within $$1.2\ensuremath{^\circ}$$',
     r'within $1.2^\circ$'),
    # Abstract and intro degree signs
    (r'$203.5\ensuremath{^\circ}$', r'$203.5^\circ$'),
    (r'$197\ensuremath{^\circ}$',   r'$197^\circ$'),
    (r'$6.5\ensuremath{^\circ}$',   r'$6.5^\circ$'),
    (r'$84.3\ensuremath{^\circ}$',  r'$84.3^\circ$'),
    (r'$92.5\ensuremath{^\circ}$',  r'$92.5^\circ$'),
    (r'$92.49\ensuremath{^\circ}$', r'$92.49^\circ$'),
    (r'$\pi/2\ensuremath{^\circ}$', r'$\pi/2$'),
    # Any remaining ensuremath degree
    (r'\ensuremath{^\circ}', r'^\circ'),
]

count = 0
for old, new in fixes:
    if old in tex:
        tex = tex.replace(old, new)
        count += 1
        print('Fixed: ' + old[:60])

# Also fix the align block issue at line 197
# The problem is a stray math mode left open from the table
# Show context around line 197
lines = tex.split('\n')
for i, line in enumerate(lines[190:205], 191):
    print(str(i) + ': ' + repr(line[:80]))

print()
print('Fixed ' + str(count) + ' patterns')
with open(path, 'w', encoding='utf-8') as f: f.write(tex)
