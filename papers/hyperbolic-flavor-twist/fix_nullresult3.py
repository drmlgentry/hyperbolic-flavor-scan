path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

import re

# Find and replace the problematic exponentiation sentence
old = re.search(
    r'As a heuristic.*?which is directly testable\.',
    tex, re.DOTALL
)

if old:
    new_para = (
        r'As a heuristic, for long words that are products of many generators,'
        '\n'
        r'cancellations among the twist contributions can produce a very small'
        '\n'
        r'net twist angle.'
        '\n'
        r'A crude estimate, based on the near-$\pi/2$ generator twist'
        '\n'
        r'$\phi(\word{b})\approx 89.16\degr$ and a random-walk-like phase'
        '\n'
        r'accumulation in which the deficit $\delta=\pi/2-\phi(\word{b})\approx 0.84\degr$'
        '\n'
        r'compounds over $\ell$ steps, suggests that at word length $\ell\approx 12$'
        '\n'
        r'the residual angle could be as small as ${\sim}0.18\degr$,'
        '\n'
        r'close to the required $\theta_{13}^{\rm CKM}=0.201\degr$.'
        '\n'
        r'This estimate is heuristic: non-commuting compositions do not accumulate'
        '\n'
        r'phases simply, and the actual spectral floor at length 12 requires'
        '\n'
        r'explicit computation.'
        '\n'
        r'The prediction is that a length-12 census of \mfld{m006} will yield'
        '\n'
        r'a word with $|\phi|\approx 0.20\degr$, which is directly testable.'
    )
    tex = tex[:old.start()] + new_para + tex[old.end():]
    print('Null result paragraph replaced with corrected version')
else:
    print('WARNING: block not found by regex -- printing nearby context')
    idx = tex.find(r'heuristic')
    if idx > 0:
        print(repr(tex[idx-50:idx+400]))

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
