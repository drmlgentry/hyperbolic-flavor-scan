import re
path = r'C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# The previous fix created patterns like $203.5^\circ$$ (double $)
# and $$92.49^\circ$$ inside already-math contexts.
# Strategy: replace ALL degree occurrences with ^\circ uniformly,
# then fix the math mode context.

# First undo the previous bad fix -- remove all $N^\circ$ patterns
# that created double-dollar issues
# Pattern 1: inside existing math: $...$X^\circ$$ -> $...X^\circ$
tex = re.sub(r'\$(\d+\.?\d*)\^\{?\\circ\}?\$\$', r'$\1^\\circ$', tex)
tex = re.sub(r'\$\$(\d+\.?\d*)\^\{?\\circ\}?\$\$', r'$\1^\\circ$', tex)

# Pattern 2: standalone degree in text: $X^\circ$ is fine outside math
# but inside math it should be X^\circ without the dollars

# Fix table entries like: & $-$168.6^\circ$$ &
# These arose from: -168.6° in a table cell already inside $...$
tex = re.sub(r'\$([+-]?\d+\.?\d*)\^\{?\\circ\}?\$\$', r'$\1^\\circ$', tex)

# Fix inline text occurrences: = $203.5^\circ$$, 
tex = re.sub(r'= \$(\d+\.?\d*)\^\{?\\circ\}?\$\$', r'= $\\1^\\circ$', tex)

# Most reliable approach: replace all degree variants with \textdegree in text
# and ^\circ in math -- but since we cant tell context easily,
# use a simpler global replacement:
# Any sequence like $X.Y^\circ$ where it appears standalone = fine
# Any sequence like $$X.Y^\circ$$ = wrong, fix to $X.Y^\circ$

# Nuclear option: replace everything back to unicode degree, 
# then do a context-aware replacement
tex = tex.replace(r'^\circ', '\u00b0')  # undo all circ replacements
tex = tex.replace(r'$\1^\\circ$', '\u00b0')  # cleanup artifacts

# Now do it properly:
# In table cells and inline text, degree sign needs \textdegree or $X^\circ$
# Replace unicode degree with $^\circ$ (works in both text and math mode
# when wrapped properly)

# Use \ensuremath{^\circ} which works in both modes
tex = tex.replace('\u00b0', r'\ensuremath{^\circ}')
print('Degrees fixed using \\ensuremath{^\\circ}')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
