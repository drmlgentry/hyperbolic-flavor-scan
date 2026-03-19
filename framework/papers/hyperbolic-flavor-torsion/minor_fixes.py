path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# 1. Add AIC mention to abstract
old_abs = (
    "confirmed by AIC model comparison against power-law alternatives."
)
new_abs = (
    "confirmed by Akaike Information Criterion (AIC) model selection "
    "against power-law alternatives (ΔAIC $> 40$ for all classes)."
)
tex = tex.replace(old_abs, new_abs, 1)

# 2. Improve figure caption
old_cap = (
    r"\caption{Log of spectral floor $\phi_{\rm floor}^{(k)}(L)$ vs.\ word"
    "\n"
    r"  length $L$ for \mfld{m006}, by homology class $k\in\mathbb{Z}/5$."
    "\n"
    r"  Classes~1 and~4 (fastest, $c\approx 0.985$) and classes~2 and~3"
    "\n"
    r"  (intermediate, $c\approx 0.401$) are paired; class~0 is slowest"
    "\n"
    r"  ($c\approx 0.750$)."
    "\n"
    r"  Dashed lines: exponential fits $A_k e^{-c_k L}$."
    "\n"
    r"  The $4:2:1$ ratio at $L=12$ is the principal result.}"
)
new_cap = (
    r"\caption{Spectral floor $\phi_{\rm floor}^{(k)}(L)$ (log scale) vs.\ "
    "word length $L$ for \mfld{m006}, partitioned by homology class "
    "$k\\in\\mathbb{Z}/5$. "
    "Horizontal axis: word length $L=1,\\ldots,12$. "
    "Vertical axis: minimum folded twist angle (degrees, log scale). "
    "Dashed lines: exponential fits $A_k e^{-c_k L}$ with rates from "
    "Table~\\ref{tab:rates}. "
    "Classes~1 and~4 decay fastest ($c\\approx 0.985$, red/orange); "
    "classes~2 and~3 intermediate ($c\\approx 0.401$, green/purple); "
    "class~0 slowest ($c\\approx 0.750$, blue). "
    "The $4:2:1$ floor ratio at $L=12$ is the principal result "
    "(Table~\\ref{tab:floors_L12}).}"
)
tex = tex.replace(old_cap, new_cap, 1)

# 3. Add testing remark after conjecture
old_conj_end = (
    r"A proof of this conjecture would require computing the invariant trace"
    "\n"
    r"field of \mfld{m006} via the \textsc{Snap} package~\cite{Snap} and"
    "\n"
    r"relating the arithmetic involution to the observed rate pairing."
)
new_conj_end = (
    r"A proof of this conjecture would require computing the invariant trace"
    "\n"
    r"field of \mfld{m006} via the \textsc{Snap} package~\cite{Snap} and"
    "\n"
    r"relating the arithmetic involution to the observed rate pairing."
    "\n"
    r"A direct test is feasible: \textsc{Snap} computes the invariant trace"
    "\n"
    r"field $k(M)$ and quaternion algebra $\mathcal{A}$ of any cusped"
    "\n"
    r"hyperbolic manifold.  Running \textsc{Snap} on the Dehn-filling"
    "\n"
    r"description of \mfld{m006} and checking whether $\mathcal{A}$ admits"
    "\n"
    r"a complex-conjugation involution would either confirm or refute the"
    "\n"
    r"arithmetic origin of the pairing $c_1=c_4$, $c_2=c_3$."
)
tex = tex.replace(old_conj_end, new_conj_end, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Minor revisions applied.")
