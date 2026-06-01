path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Fix 1: add corollary to newtheorem block
old = r'\newtheorem{remark}{Remark}'
new = r'\newtheorem{remark}{Remark}' + '\n' + r'\newtheorem{corollary}{Corollary}'
if 'newtheorem{corollary}' not in tex:
    tex = tex.replace(old, new)
    print('Fix 1: corollary newtheorem added')
else:
    print('Fix 1: corollary already defined')

# Fix 2: replace degree symbols in table cells
# The degree symbol as unicode char causes issues in math mode
# Replace 92.487° -> $92.487^\circ$ etc.
import re
# Pattern: number followed by degree symbol (unicode °)
tex = re.sub(r'(\d+\.\d+)\u00b0', r'$\1^\\circ$', tex)
tex = re.sub(r'(\d+)\u00b0', r'$\1^\\circ$', tex)
print('Fix 2: degree symbols replaced')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
