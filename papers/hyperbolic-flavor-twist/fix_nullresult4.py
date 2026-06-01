path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

import re

old = re.search(
    r'\\section\{Null Result:.*?\\label\{sec:null\}.*?(?=\\section)',
    tex, re.DOTALL
)

if old:
    new_section = r"""\section{Null Result: \texorpdfstring{$\theta_{13}^{\rm CKM}$}{theta13(CKM)}}
\label{sec:null}

Six of seven independent SM flavor parameters are matched to better than $3\%$.
The sole outlier is $\theta_{13}^{\rm CKM}=0.201\degr$.
The minimum $|\phi|$ at length~5 on \mfld{m006} is
$\phi(\word{abbAb})=1.687\degr$, a discrepancy of $1.49\degr$ ($740\%$).

The angle $\theta_{13}^{\rm CKM}$ is parametrically of order
$\lambda_W^3\approx 0.011$ in the Wolfenstein parameterisation,
making it the smallest and most suppressed of the quark mixing angles.

For long words composed of many generators, cancellations among the
twist contributions can produce a very small net twist angle; a
crude estimate based on the near-$\pi/2$ generator twist
$\phi(\word{b})\approx 89.16\degr$ and the deficit
$\delta=\pi/2-\phi(\word{b})\approx 0.84\degr$ compounding over
$\ell$ steps suggests a spectral floor of ${\sim}0.18\degr$ at
$\ell\approx 12$.
This estimate is heuristic: non-commuting compositions do not accumulate
phases simply, and the actual spectral floor at length 12 requires
explicit computation.

A rough extrapolation places $\theta_{13}^{\rm CKM}$ within reach at
$\ell^*\approx 12$--$14$, constituting a falsifiable prediction of the
framework: a length-12 census of \mfld{m006} should yield a word with
$|\phi|\approx 0.20\degr$.

This early appearance of SM parameters---all but one by word length
five---is consistent with the large spectral gap $\lambda_1\sim 2.5$--$2.8$
established in Section~\ref{sec:selberg}: the gap forces the low-length
geodesic spectrum to be sufficiently rich and well-distributed that the SM
angles emerge at the shortest scales, rather than being buried deep in the
length spectrum where coincidences would be less significant.

"""
    tex = tex[:old.start()] + new_section + tex[old.end():]
    print('Null result section replaced')
else:
    print('WARNING: section not found by regex -- printing context')
    idx = tex.find('sec:null')
    if idx > 0:
        print(repr(tex[idx-50:idx+500]))

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
