path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Fix 1: suppress hyperref math warnings with texorpdfstring in remaining headings
tex = tex.replace(
    r'\section{Null Result: $\theta_{13}^{\rm CKM}$}',
    r'\section{Null Result: \texorpdfstring{$\theta_{13}^{\rm CKM}$}{theta13(CKM)}}'
)
tex = tex.replace(
    r'\section{The $M_Z/M_W$ Ratio}',
    r'\section{The \texorpdfstring{$M_Z/M_W$}{MZ/MW} Ratio}'
)
tex = tex.replace(
    r'\section{The $\overline{\rm MS}$ quark mass ratio}',
    r'\section{The \texorpdfstring{$\overline{\mathrm{MS}}$}{MSbar} quark mass ratio}'
)
tex = tex.replace(
    r'\subsection{Numerical estimates of $\lambda_1$}',
    r'\subsection{Numerical estimates of \texorpdfstring{$\lambda_1$}{lambda1}}'
)
print('Fix 1: texorpdfstring restored with safe strings')

# Fix 2: fix overfull longtable row -- shorten the caption and column widths
# The issue is the longtable header line is too wide
# Add @{} to remove inter-column padding and tighten column spec
tex = tex.replace(
    r'\begin{longtable}{lllllr}',
    r'\begin{longtable}{@{}lllllr@{}}'
)
print('Fix 2: longtable column padding tightened')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
