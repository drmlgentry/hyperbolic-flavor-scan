path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Fix 1: remove \orcid line -- not supported in all elsarticle versions
tex = tex.replace(
    r'\orcid{0000-0009-0006-4550-2663}' + '\n',
    ''
)
print('Fix 1: orcid removed')

# Fix 2: switch from authoryear to numbered citations
tex = tex.replace(
    r'\documentclass[preprint,12pt,authoryear]{elsarticle}',
    r'\documentclass[preprint,12pt,number]{elsarticle}'
)
print('Fix 2: citation style -> number')

# Fix 3: \journal must be INSIDE \begin{document}, not before it
# Move it to after \begin{document}
tex = tex.replace(
    r'\journal{Nuclear Physics B}' + '\n\n' + r'\begin{document}',
    r'\begin{document}' + '\n' + r'\journal{Nuclear Physics B}' + '\n'
)
print('Fix 3: journal macro moved inside document')

# Fix 4: add ORCID as a note instead
tex = tex.replace(
    r'\ead{drmlgentry@protonmail.com}' + '\n' + r'\address[a]',
    r'\ead{drmlgentry@protonmail.com}' + '\n'
    r'\ead[orcid]{https://orcid.org/0009-0006-4550-2663}' + '\n'
    r'\address[a]'
)
print('Fix 4: ORCID as ead url')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
