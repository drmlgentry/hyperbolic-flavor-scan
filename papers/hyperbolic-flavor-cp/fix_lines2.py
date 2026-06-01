path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: lines = f.readlines()

# Fix line 158: PMNS table header (currently has escaped backslashes)
lines[157] = 'Word & $t$ & $\\phi$ (deg) & $|\\lambda| = e^{|t|}$ \\\\\n'

# Fix line 486: CKM table header (corrupted modulus text)
lines[485] = 'Word & $t$ & $\\phi$ (deg) & $|\\lambda| = e^{|t|}$ \\\\\n'

print('Line 158 now: ' + repr(lines[157]))
print('Line 486 now: ' + repr(lines[485]))

with open(path, 'w', encoding='utf-8') as f: f.writelines(lines)
print('Done')
