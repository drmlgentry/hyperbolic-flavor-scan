path_in  = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist.tex'
path_out = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'

with open(path_in, 'r', encoding='utf-8-sig') as f: tex = f.read()

# 1. Replace documentclass + jheppub with elsarticle
tex = tex.replace(
    r'\documentclass[a4paper,11pt]{article}',
    r'\documentclass[preprint,12pt,authoryear]{elsarticle}'
)
tex = tex.replace(
    r'\usepackage{jheppub}',
    r'\usepackage{lineno}'
)

# 2. Replace title block (jhep style -> elsarticle style)
old_title_block = (
    r'\title{\boldmath Twist Angle Spectrum of Hyperbolic Holonomy:\\' + '\n'
    r'Encoding of Standard Model Flavor Parameters}' + '\n\n'
    r'\author{Marvin L.~Gentry}' + '\n'
    r'\affiliation{Independent Researcher, Seattle, WA, USA}' + '\n'
    r'\emailAdd{drmlgentry@protonmail.com}' + '\n\n'
    r'\abstract{'
)
new_title_block = (
    r'\begin{frontmatter}' + '\n\n'
    r'\title{Twist Angle Spectrum of Hyperbolic Holonomy:\\' + '\n'
    r'Encoding of Standard Model Flavor Parameters}' + '\n\n'
    r'\author[a]{Marvin L.~Gentry}' + '\n'
    r'\ead{drmlgentry@protonmail.com}' + '\n'
    r'\address[a]{Independent Researcher, Seattle, WA, USA}' + '\n\n'
    r'\begin{abstract}'
)
if old_title_block in tex:
    tex = tex.replace(old_title_block, new_title_block, 1)
    print('Title block converted')
else:
    print('WARNING: title block not found -- check manually')

# 3. Close abstract and add keywords, close frontmatter
old_abstract_end = (
    r'\keywords{hyperbolic 3-manifolds, holonomy, loxodromic elements,' + '\n'
    r'Standard Model flavor, CKM, PMNS, quark masses, Selberg zeta function,' + '\n'
    r'Iwasawa decomposition}'
)
new_abstract_end = (
    r'\end{abstract}' + '\n\n'
    r'\begin{keyword}' + '\n'
    r'hyperbolic 3-manifolds \sep holonomy \sep loxodromic elements \sep' + '\n'
    r'Standard Model flavor \sep CKM \sep PMNS \sep quark masses \sep' + '\n'
    r'Selberg zeta function \sep Iwasawa decomposition' + '\n'
    r'\end{keyword}' + '\n\n'
    r'\end{frontmatter}' + '\n\n'
    r'\linenumbers'
)
if old_abstract_end in tex:
    tex = tex.replace(old_abstract_end, new_abstract_end, 1)
    print('Abstract/keywords converted')
else:
    print('WARNING: keywords block not found -- check manually')

# 4. Replace \begin{document}\maketitle with just \begin{document}
tex = tex.replace(
    r'\begin{document}' + '\n' + r'\maketitle',
    r'\begin{document}'
)

# 5. Replace acknowledgments environment
tex = tex.replace(r'\acknowledgments', r'\section*{Acknowledgements}')

# 6. Replace thebibliography with natbib-style
# elsarticle uses natbib; keep thebibliography as-is but switch \bibitem style
# Just need to ensure \bibliographystyle is not needed (inline bib is fine)
# Add orcid info after email
tex = tex.replace(
    r'\ead{drmlgentry@protonmail.com}',
    r'\ead{drmlgentry@protonmail.com}' + '\n'
    r'\orcid{0000-0009-0006-4550-2663}'
)

# 7. Fix \cite style -- elsarticle/natbib uses \cite{} fine as-is
# but \citeauthor etc not used here so nothing to fix

# 8. Add journal macro for NPB
npb_macro = (
    r'\journal{Nuclear Physics B}' + '\n\n'
)
tex = tex.replace(
    r'\begin{document}',
    npb_macro + r'\begin{document}',
    1
)
print('Journal macro added')

with open(path_out, 'w', encoding='utf-8') as f: f.write(tex)
print('NPB version written to: ' + path_out)
