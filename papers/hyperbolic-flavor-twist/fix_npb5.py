path = r'C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-npb.tex'
with open(path, 'r', encoding='utf-8-sig') as f: tex = f.read()

# Fix 1: remove the broken \small{ group wrapper
tex = tex.replace('{\small\n\\begin{longtable}', '\\begin{longtable}')
tex = tex.replace('\\end{longtable}\n}', '\\end{longtable}', 1)
print('Fix 1: small group removed')

# Fix 2: replace the entire longtable block with a compact regular table
import re
old = re.search(
    r'\\begin\{longtable\}.*?\\end\{longtable\}',
    tex, re.DOTALL
)
if old:
    new_table = r"""\begin{table}[htbp]
\centering
\caption{Complete sub-$5\degr$ and sub-$2\%$ A-factor spectral coincidences
  up to word length~5. Fold: d$=$direct, $\pi$$=180\degr-|\phi|$,
  r$=$modulus ratio. Length$=\max(|W_1|,|W_2|)$ for ratio rows.}
\label{tab:full}
\small
\begin{tabular}{@{}llp{1.6cm}p{1.6cm}lr@{}}
\toprule
Mfld & Word & $|\phi|$ (deg) & Fold & SM match & Error \\
\midrule
\multicolumn{6}{@{}l}{\textit{$\delta_{\rm CKM}=68.0\degr$}}\\
\mfld{m003} & \word{bbbb} & $67.428$ & $\pi$ & $\delta_{\rm CKM}$ & $0.57\degr$ \\
\mfld{m006} & \word{aa} & $67.654$ & $\pi$ & $\delta_{\rm CKM}$ & $0.35\degr$ \\
\midrule
\multicolumn{6}{@{}l}{\textit{$\theta_{12}^\nu=33.41\degr$}}\\
\mfld{m003} & \word{ABaB} & $32.980$ & d & $\theta_{12}^\nu$ & $0.43\degr$ \\
\mfld{m006} & \word{abbAB} & $33.620$ & $\pi$ & $\theta_{12}^\nu$ & $0.21\degr$ \\
\midrule
\multicolumn{6}{@{}l}{\textit{$\theta_{23}^\nu=49.1\degr$}}\\
\mfld{m003} & \word{aabb} & $48.081$ & $\pi$ & $\theta_{23}^\nu$ & $1.02\degr$ \\
\midrule
\multicolumn{6}{@{}l}{\textit{$\theta_{12}^{\rm CKM}=13.04\degr$}}\\
\mfld{m003} & \word{AAB} & $12.638$ & $\pi$ & $\theta_{12}^{\rm CKM}$ & $0.40\degr$ \\
\midrule
\multicolumn{6}{@{}l}{\textit{$7\theta_{12}^{\rm CKM}=91.28\degr$}}\\
\mfld{m003} & \word{AAAB} & $90.750$ & d & $7\theta_{12}^{\rm CKM}$ & $0.53\degr$ \\
\mfld{m006} & \word{AAb} & $92.487$ & d & $7\theta_{12}^{\rm CKM}$ & $1.21\degr$ \\
\midrule
\multicolumn{6}{@{}l}{\textit{$\bar{m}_b/\bar{m}_c=3.2913$ (modulus ratios)}}\\
\mfld{m003} & \word{bbbb}/\word{bAbA} & $3.2910$ & r & $\bar{m}_b/\bar{m}_c$ & $0.01\%$ \\
\mfld{m003} & \word{bbbb}/\word{AA} & $3.2910$ & r & $\bar{m}_b/\bar{m}_c$ & $0.01\%$ \\
\midrule
\multicolumn{6}{@{}l}{\textit{$M_Z/M_W=1.13451$ (modulus ratios)}}\\
\mfld{m006} & \word{aaabb}/\word{baaab} & $1.1355$ & r & $M_Z/M_W$ & $0.09\%$ \\
\bottomrule
\end{tabular}
\end{table}"""
    tex = tex[:old.start()] + new_table + tex[old.end():]
    print('Fix 2: longtable replaced with compact table')
else:
    print('WARNING: longtable not found')

with open(path, 'w', encoding='utf-8') as f: f.write(tex)
print('Done')
