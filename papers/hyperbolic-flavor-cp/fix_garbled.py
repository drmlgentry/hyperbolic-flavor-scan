path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: lines = f.readlines()

# Show all corrupted lines first
print("Scanning for corrupted lines...")
for i, line in enumerate(lines, 1):
    if 'modulust' in line or 'modulusi' in line or 'modulus|t' in line or 'modulus^' in line:
        print(str(i) + ': ' + repr(line))

# Fix each one
fixed = 0
for i, line in enumerate(lines):
    # Introduction corruption: "modulus" concatenated with "t(\gamma)}"
    if 'modulust(\\gamma)}' in line or 'modulus^{|t(\\gamma)|}' in line:
        lines[i] = '$e^{|t(\\gamma)|}$, while the full complex eigenvalue\n'
        print('Fixed intro line ' + str(i+1))
        fixed += 1

    # Discussion corruption: "modulusi\phi}"
    if 'modulusi\\phi}' in line or 'moduluse^{i' in line:
        lines[i] = 'is notable: $e^{i\\phi} \\approx i$, so the holonomy eigenvalue\n'
        print('Fixed discussion phi line ' + str(i+1))
        fixed += 1

    # Discussion corruption: "modulus|t|}"
    if 'modulus|t|}' in line:
        lines[i] = 'The translation lengths $e^{|t|}$ of the three optimal PMNS words\n'
        print('Fixed discussion |t| line ' + str(i+1))
        fixed += 1

    # Conclusion corruption: same as introduction
    if 'modulust(\\gamma' in line and i > 400:
        lines[i] = '$e^{|t(\\gamma)|}$, while the full complex eigenvalue\n'
        print('Fixed conclusion line ' + str(i+1))
        fixed += 1

print(str(fixed) + ' lines fixed')
with open(path, 'w', encoding='utf-8') as f: f.writelines(lines)
