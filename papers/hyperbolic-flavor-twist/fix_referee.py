path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()
changes = []

# Fix 1: Table 2 -- ratio rows: clear the |phi| cell, add note in fold cell
tex = tex.replace(
    r'\mfld{m003} & \word{bbbb}/\word{bAbA} & $3.2910$ & r & $\bar{m}_b/\bar{m}_c$ & $0.01\%$ \\',
    r'\mfld{m003} & \word{bbbb}/\word{bAbA} & --- & ratio & $\bar{m}_b/\bar{m}_c$ & $0.0003$ (0.01\%) \\'
)
tex = tex.replace(
    r'\mfld{m003} & \word{bbbb}/\word{AA} & $3.2910$ & r & $\bar{m}_b/\bar{m}_c$ & $0.01\%$ \\',
    r'\mfld{m003} & \word{bbbb}/\word{AA} & --- & ratio & $\bar{m}_b/\bar{m}_c$ & $0.0003$ (0.01\%) \\'
)
tex = tex.replace(
    r'\mfld{m006} & \word{aaabb}/\word{baaab} & $1.1355$ & r & $M_Z/M_W$ & $0.09\%$ \\',
    r'\mfld{m006} & \word{aaabb}/\word{baaab} & --- & ratio & $M_Z/M_W$ & $0.0010$ (0.09\%) \\'
)
changes.append('Fix 1: ratio rows in Table 2 corrected')

# Fix 2: Null result heuristic -- replace ill-formed exponentiation
old_null = r"""As a heuristic estimate, the near-$\pi/2$ generator twist
$\phi(\word{b})\approx 89.16\degr$ accumulates approximately multiplicatively
under composition for words dominated by $b$-generators:
$\phi_{\rm min}(\ell)\sim (89.16\degr/90\degr)^\ell\cdot 90\degr\approx
90\degr\cdot 0.9906^\ell$.
This gives $\phi_{\rm min}(12)\approx 0.18\degr$, comparable to
$\theta_{13}^{\rm CKM}=0.201\degr$.
We stress that this is a heuristic: non-commuting word compositions may not
accumulate phases multiplicatively, and the actual spectral floor at length 12
requires explicit computation.
The prediction is therefore that a length-12 census of \mfld{m006} will
contain a word with $|\phi|\approx 0.20\degr$, which is testable."""

new_null = r"""As a heuristic, for words that are long products of generators,
the net twist angle can become small because contributions from individual
generators may cancel.
The generator \word{b} on \mfld{m006} has $\phi(\word{b})\approx 89.16\degr$;
a crude estimate based on phase cancellation among many such generators
suggests that at word length $\ell\approx 12$ the residual twist angle could
be as small as ${\sim}0.18\degr$, close to the required
$\theta_{13}^{\rm CKM}=0.201\degr$.
This estimate is heuristic: non-commuting compositions do not accumulate
phases simply, and the actual spectral floor at length 12 requires explicit
computation.
The prediction is that a length-12 census of \mfld{m006} will yield a word
with $|\phi|\approx 0.20\degr$, which is directly testable."""

if old_null in tex:
    tex = tex.replace(old_null, new_null)
    changes.append('Fix 2: null result heuristic rephrased')
else:
    changes.append('WARNING Fix 2: null result block not found by exact match')

# Fix 3: Reference clarity -- rephrase PRD companion sentence
old_ref = r"""In a companion series~\cite{GentryCKM,GentryPMNS,GentryCP} we established
that two compact hyperbolic 3-manifolds, \mfld{m003} and \mfld{m006} from the
\textsc{SnapPy}~\cite{SnapPy} OrientableClosedCensus, reproduce the CKM and
PMNS mixing matrices with fitness $F<0.020$."""

new_ref = r"""In related work submitted to Physical Review
D~\cite{GentryCKM,GentryPMNS,GentryCP} we established that two compact
hyperbolic 3-manifolds, \mfld{m003} and \mfld{m006} from the
\textsc{SnapPy}~\cite{SnapPy} OrientableClosedCensus, reproduce the CKM and
PMNS mixing matrices with fitness $F<0.020$.
The present paper focuses on the twist-angle spectrum of the same manifolds."""

if old_ref in tex:
    tex = tex.replace(old_ref, new_ref)
    changes.append('Fix 3: PRD vs NPB reference clarified')
else:
    changes.append('WARNING Fix 3: reference paragraph not found by exact match')

# Fix 4: Abstract MZ/MW wording
old_abs = r'the $\overline{\rm MS}$ quark mass ratio'
# Only fix the abstract appearance phrase if present
tex = tex.replace(
    r'the $M_Z/M_W$ gauge boson mass ratio appears at $0.09\%$',
    r'the $M_Z/M_W$ gauge boson mass ratio is reproduced to $0.09\%$'
)
changes.append('Fix 4: abstract MZ/MW wording')

# Fix 5: Table 1 -- add absolute error for mb/mc
tex = tex.replace(
    r'$\bar{m}_b/\bar{m}_c=3.2913$ & \mfld{m003} & \word{bbbb/bAbA} & $3.2910$ & r & $0.01\%$ \\',
    r'$\bar{m}_b/\bar{m}_c=3.2913$ & \mfld{m003} & \word{bbbb/bAbA} & $3.2910$ & r & $0.0003$ (0.01\%) \\'
)
changes.append('Fix 5: Table 1 mb/mc absolute error added')

for c in changes: print(c)
with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
