path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

old = r"""\small
\begin{tabular}{llllllr}
\toprule
SM Parameter & Target & Mfld & Word & $|\phi|$ or ratio & Fold & Error \\
\midrule
$\theta_{12}^{\rm CKM}=13.04\degr$ & $13.04\degr$ & \mfld{m003} & \word{AAB} & $167.362\degr$ & $180\degr$ & $0.402\degr$ (3.1\%) \\
$\theta_{23}^{\rm CKM}=2.38\degr$ & $2.38\degr$ & \mfld{m006} & \word{aaabb} & $2.132\degr$ & direct & $0.248\degr$ (10.4\%) \\
$\delta_{\rm CKM}=68.0\degr$ & $68.0\degr$ & \mfld{m006} & \word{aa} & $112.346\degr$ & $180\degr$ & $0.346\degr$ (0.51\%) \\
$\theta_{12}^\nu=33.41\degr$ & $33.41\degr$ & \mfld{m006} & \word{abbAB} & $146.380\degr$ & $180\degr$ & $0.210\degr$ (0.63\%) \\
$\theta_{13}^\nu=8.54\degr$ & $8.54\degr$ & \mfld{m003} & \word{AA} & $168.556\degr$ & $180\degr$ & $2.904\degr$ (34\%) \\
$\theta_{23}^\nu=49.1\degr$ & $49.1\degr$ & \mfld{m003} & \word{aabb} & $131.919\degr$ & $180\degr$ & $1.019\degr$ (2.1\%) \\
$\delta_{\rm CP}=197\degr$ & $197\degr$ & \mfld{m003} & \word{aa/ab/aB} & signed sum & --- & $6.5\degr$ (3.3\%) \\
$m_b/m_c=3.291$ & $3.2913$ & \mfld{m003} & \word{bbbb/bAbA} & $3.2910$ & ratio & $0.0003$ (0.01\%) \\
$M_Z/M_W=1.1345$ & $1.1345$ & \mfld{m006} & \word{aaabb/baaab} & $1.1355$ & ratio & $0.0010$ (0.09\%) \\
\bottomrule
\end{tabular}"""

new = r"""\small
\begin{tabular}{@{}lllp{1.5cm}lr@{}}
\toprule
SM Parameter & Mfld & Word & $|\phi|$/ratio & Fold & Error \\
\midrule
$\theta_{12}^{\rm CKM}=13.04\degr$ & \mfld{m003} & \word{AAB} & $167.362\degr$ & $180\degr$ & $0.40\degr$ (3.1\%) \\
$\theta_{23}^{\rm CKM}=2.38\degr$ & \mfld{m006} & \word{aaabb} & $2.132\degr$ & d & $0.25\degr$ (10.4\%) \\
$\delta_{\rm CKM}=68.0\degr$ & \mfld{m006} & \word{aa} & $112.346\degr$ & $180\degr$ & $0.35\degr$ (0.51\%) \\
$\theta_{12}^\nu=33.41\degr$ & \mfld{m006} & \word{abbAB} & $146.380\degr$ & $180\degr$ & $0.21\degr$ (0.63\%) \\
$\theta_{13}^\nu=8.54\degr$ & \mfld{m003} & \word{AA} & $168.556\degr$ & $180\degr$ & $2.90\degr$ (34\%) \\
$\theta_{23}^\nu=49.1\degr$ & \mfld{m003} & \word{aabb} & $131.919\degr$ & $180\degr$ & $1.02\degr$ (2.1\%) \\
$\delta_{\rm CP}=197\degr$ & \mfld{m003} & \word{aa/ab/aB} & signed sum & --- & $6.5\degr$ (3.3\%) \\
$\bar{m}_b/\bar{m}_c=3.2913$ & \mfld{m003} & \word{bbbb/bAbA} & $3.2910$ & r & $0.01\%$ \\
$M_Z/M_W=1.13451$ & \mfld{m006} & \word{aaabb/baaab} & $1.1355$ & r & $0.09\%$ \\
\bottomrule
\end{tabular}"""

if old in tex:
    tex = tex.replace(old, new)
    print('Table replaced')
else:
    print('WARNING: old table not found -- check whitespace')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
