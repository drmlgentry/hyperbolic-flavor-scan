path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Append the spectral gap paragraph after the existing testable prediction sentence
old_end = r'with $|\phi|\approx 0.20\degr$, which is directly testable.'

new_end = (
    r'with $|\phi|\approx 0.20\degr$, which is directly testable.' + '\n\n'
    r'This early appearance of SM parameters---all but one by word length' + '\n'
    r'five---is consistent with the large spectral gap $\lambda_1\sim 2.5$--$2.8$' + '\n'
    r'established in Section~\ref{sec:selberg}: the gap forces the low-length' + '\n'
    r'geodesic spectrum to be sufficiently rich and well-distributed that the SM' + '\n'
    r'angles emerge at the shortest scales, rather than being buried deep in the' + '\n'
    r'length spectrum where coincidences would be less significant.'
)

if old_end in tex:
    tex = tex.replace(old_end, new_end, 1)
    print('Spectral gap paragraph appended after null result')
else:
    print('WARNING: anchor sentence not found -- paste surrounding text')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
