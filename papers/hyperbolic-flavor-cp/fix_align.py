path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: lines = f.readlines()

fixes = {
    # Line 194: inside align, no $ needed around numbers
    193: r'\delta_{\mathrm{geom}} &= (-168.6^\circ) - (83.7^\circ) + (95.7^\circ) \nonumber \\' + '\n',
    # Line 195: inside align
    194: r'&= -156.6^\circ \equiv 203.5^\circ \pmod{360^\circ}.' + '\n',
    # Line 198: inline text, degrees in math mode
    197: r'The PDG 2024 central value is $\delta_{\mathrm{CP}} = 197^\circ$,' + '\n',
    # Line 199: stray $$ 
    198: r'giving a discrepancy of $6.5^\circ$ or $3.3\%$.' + '\n',
}

for idx, new_line in fixes.items():
    print('Line ' + str(idx+1) + ' was: ' + repr(lines[idx].strip()))
    lines[idx] = new_line
    print('Line ' + str(idx+1) + ' now: ' + repr(lines[idx].strip()))
    print()

with open(path, 'w', encoding='utf-8') as f: f.writelines(lines)
print('Done')
