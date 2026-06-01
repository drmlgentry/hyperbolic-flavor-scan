path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Fix 1: add hyperref (provides \texorpdfstring and \url)
tex = tex.replace(
    r'\usepackage{lineno}',
    r'\usepackage{lineno}' + '\n' + r'\usepackage{hyperref}'
)
print('Fix 1: hyperref added')

# Fix 2: replace \texorpdfstring in section headings with plain text
tex = tex.replace(
    r'\section{The \texorpdfstring{$\overline{\rm MS}$}{MSbar} quark mass ratio}',
    r'\section{The $\overline{\rm MS}$ quark mass ratio}'
)
tex = tex.replace(
    r'\subsection{The \texorpdfstring{$\overline{\rm MS}$}{MSbar} quark mass ratio}',
    r'\subsection{The $\overline{\rm MS}$ quark mass ratio}'
)
tex = tex.replace(
    r'\section{The \texorpdfstring{$M_Z/M_W$}{MZ/MW} Ratio}',
    r'\section{The $M_Z/M_W$ Ratio}'
)
tex = tex.replace(
    r'\section{Null Result: \texorpdfstring{$\theta_{13}^{\rm CKM}$}{theta13 CKM}}',
    r'\section{Null Result: $\theta_{13}^{\rm CKM}$}'
)
tex = tex.replace(
    r'\subsection{Numerical estimates of \texorpdfstring{$\lambda_1$}{lambda1}}',
    r'\subsection{Numerical estimates of $\lambda_1$}'
)
tex = tex.replace(
    r'\subsection{Large spectral gap and SM-scale mixing}',
    r'\subsection{Large spectral gap and SM-scale mixing}'
)
print('Fix 2: texorpdfstring replaced')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
