import re

path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    content = f.read()

# Find non-ASCII characters and report their positions
bad = []
for i, ch in enumerate(content):
    if ord(ch) > 127:
        line_num = content[:i].count('\n') + 1
        bad.append((line_num, i, repr(ch), ord(ch)))

if bad:
    print(f"Found {len(bad)} non-ASCII characters:")
    for line_num, pos, ch, code in bad[:20]:
        context = content[max(0,pos-20):pos+20].replace('\n','|')
        print(f"  Line {line_num}: U+{code:04X} {ch}  ...{context}...")
else:
    print("No non-ASCII characters found.")
