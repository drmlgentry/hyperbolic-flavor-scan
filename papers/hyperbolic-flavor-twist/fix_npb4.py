path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# The overfull is in the longtable caption line -- it's too long for preprint width
# Replace longtable with p{} columns to force wrapping
tex = tex.replace(
    r'\begin{longtable}{@{}lllllr@{}}',
    r'\begin{longtable}{@{}p{0.8cm}p{2.2cm}p{0.4cm}p{1.8cm}p{0.8cm}r@{}}'
)

# Also shorten the endfirsthead/endhead header lines which repeat
# The caption itself may be too wide -- add \small before longtable
tex = tex.replace(
    r'\begin{longtable}{@{}p{0.8cm}',
    r'{\small' + '\n' + r'\begin{longtable}{@{}p{0.8cm}'
)
# Close the \small group after \end{longtable}
tex = tex.replace(
    r'\end{longtable}',
    r'\end{longtable}' + '\n' + r'}'
, 1)  # only first occurrence

print('Fix: longtable column widths set explicitly')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
