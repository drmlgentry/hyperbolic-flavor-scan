path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: lines = f.readlines()

# Show current state of problem lines
for i in [157, 158, 484, 485, 486, 487, 488]:
    print(str(i+1) + ': ' + repr(lines[i]))

# Fix line 158: table header with escaped backslashes
lines[157] = 'Word & $t$ & $\\phi$ (deg) & $|\\lambda| = e^{|t|}$ \\\\\n'

# Find and fix the second broken table header (CKM table around line 487)
for i, line in enumerate(lines):
    if 'midrule' in line and i > 480 and i < 500:
        print('midrule at line ' + str(i+1))
        print('  prev: ' + repr(lines[i-1]))
        
with open(path, 'w', encoding='utf-8') as f: f.writelines(lines)
print('Fixed line 158')
