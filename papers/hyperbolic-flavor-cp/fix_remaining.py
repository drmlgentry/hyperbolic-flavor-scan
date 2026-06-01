path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

import re

# Find ALL remaining broken degree patterns and show them
print("Scanning for broken degree patterns...")
lines = tex.split('\n')
for i, line in enumerate(lines, 1):
    if '^\circ' in line or '\u00b0' in line:
        print(str(i) + ': ' + repr(line[:100]))

print()
print("Applying fixes...")

# Fix every pattern of the form:  $X^\circ$  where the $ are mismatched
# The safe approach: replace ALL degree occurrences with a single consistent form
# Rule: if already inside $...$ context, use ^\circ
# If in text, use $^\circ$
# 
# Since we can't reliably detect context, use this strategy:
# 1. Replace every  $NUMBER^\circ$  ->  NUMBER^\circ  (strip outer dollars)
# 2. Then re-wrap each number+circ in single $...$
# 3. Special: $-NUMBER^\circ$ and $+NUMBER^\circ$ patterns

# Step 1: normalize - strip all dollar signs around degree values
# Pattern: optional $ then optional sign then digits then ^\circ then optional $
tex = re.sub(r'\$([+-]?)(\d+\.?\d*)\^\\circ\$', r'\1\2^\\circ', tex)
# Double dollar variants
tex = re.sub(r'\$\$([+-]?)(\d+\.?\d*)\^\\circ\$\$', r'\1\2^\\circ', tex)
tex = re.sub(r'\$\$([+-]?)(\d+\.?\d*)\^\\circ\$', r'\1\2^\\circ', tex)
tex = re.sub(r'\$([+-]?)(\d+\.?\d*)\^\\circ\$\$', r'\1\2^\\circ', tex)

print("Pass 1 done (strip dollars)")

# Step 2: now wrap each bare NUMBER^\circ in $...$
# But NOT if already inside a math environment (align, equation)
# Simple heuristic: wrap all of them, align/equation won't mind extra $ if already in math
# Actually in align we should NOT add $, in text we should
# Use \ensuremath which works in both:
tex = re.sub(r'([+-]?\d+\.?\d*)\^\\circ', r'\\ensuremath{\1^\\circ}', tex)
print("Pass 2 done (wrap with ensuremath)")

# Step 3: fix specific known bad patterns that remain
remaining = [
    # abstract line 43
    (r'= \ensuremath{203.5^\circ},',      r'= $203.5^\circ$,'),
    (r'of \ensuremath{197^\circ} with',   r'of $197^\circ$ with'),
    # intro line 93
    (r'= \ensuremath{203.5^\circ}$',      r'= 203.5^\circ$'),
    # phi_a line 213  
    (r'= -\ensuremath{84.3^\circ}$',      r'= -84.3^\circ$'),
    # CKM section line 461
    (r'$7 \times \ensuremath{13.04^\circ} = \ensuremath{91.28^\circ}$,',
     r'$7 \times 13.04^\circ = 91.28^\circ$,'),
    (r'within \ensuremath{1.2^\circ}$',   r'within $1.2^\circ$'),
    (r'within $\ensuremath{1.2^\circ}$',  r'within $1.2^\circ$'),
    # conclusion line 486
    (r'= \ensuremath{203.5^\circ}$,',     r'= 203.5^\circ$,'),
    (r'of \ensuremath{197^\circ}$,',      r'of $197^\circ$,'),
    # pmod in align
    (r'\pmod{\ensuremath{360^\circ}}',    r'\pmod{360^\circ}'),
]
for old, new in remaining:
    if old in tex:
        tex = tex.replace(old, new)
        print('Fixed: ' + old[:60])

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print()
print("Done. Remaining ^\circ occurrences:")
for i, line in enumerate(tex.split('\n'), 1):
    if '^\circ' in line:
        print(str(i) + ': ' + repr(line[:100]))
