path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# ── Fix 1: Reframe introduction ────────────────────────────────────
old_intro = (
    r"In a companion paper~\cite{GentryTwist} we performed a systematic census of"
    "\n"
    r"loxodromic twist angles $\phi(\gamma)=\mathrm{Im}\log\lambda(\gamma)$ for all"
    "\n"
    r"closed geodesics up to word length~7 in \mfld{m003} and \mfld{m006}, two"
    "\n"
    r"compact hyperbolic 3-manifolds that reproduce the Standard Model (SM) CKM and"
    "\n"
    r"PMNS mixing matrices~\cite{GentryCKM,GentryPMNS}."
    "\n"
    r"A striking feature emerged: the A-factor spectrum of \mfld{m006} contains a"
    "\n"
    r"clean gap $[0.005\degr, 1.61\degr]$ at lengths $\leq 7$, within which the"
    "\n"
    r"CKM reactor angle $\theta_{13}^{\rm CKM}=0.201\degr$ lies."
    "\n\n"
    r"The origin of this gap was left unexplained."
    "\n"
    r"Here we show that it arises from the \emph{torsion structure} of"
    "\n"
    r"$H_1(\mfld{m006},\mathbb{Z})=\mathbb{Z}/5$."
    "\n"
    r"The fundamental group $\pi_1(M)$ maps onto $H_1(M,\mathbb{Z})$ by abelianization,"
    "\n"
    r"partitioning closed geodesics into 5 cosets."
    "\n"
    r"We compute the spectral floor --- the minimum $|\phi_{\rm fold}|$ --- for each"
    "\n"
    r"coset and find that coset~4 on \mfld{m006} has a floor 637 times smaller than"
    "\n"
    r"the next lowest coset."
    "\n"
    r"The near-identity geodesic $\word{AAABAB}$ (length~6, $\phi_{\rm fold}=0.0053\degr$)"
    "\n"
    r"belongs exclusively to coset~4, providing a concrete example of a"
    "\n"
    r"\emph{torsion-obstructed} geodesic: non-contractible (non-trivial in $H_1$) yet"
    "\n"
    r"with holonomy nearly equal to the identity."
)

new_intro = (
    r"The complex length spectrum $\{{\ell(\gamma)+i\phi(\gamma)}\}_\gamma$ of a"
    "\n"
    r"compact hyperbolic 3-manifold $M$ is a fundamental geometric invariant,"
    "\n"
    r"encoding both the geodesic lengths $\ell(\gamma)$ and loxodromic twist"
    "\n"
    r"angles $\phi(\gamma)=\mathrm{Im}\log\lambda(\gamma)$."
    "\n"
    r"While the length spectrum is well-studied~\cite{MaclachlanReid},"
    "\n"
    r"the distribution of twist angles and their relationship to the homological"
    "\n"
    r"structure of $M$ has received less attention."
    "\n\n"
    r"In this paper we study the partition of the twist angle spectrum"
    "\n"
    r"$\Phi(M)=\{|\phi_{\rm fold}(\gamma)|\}_\gamma$ by the abelianization map"
    "\n"
    r"$\pi_1(M)\to H_1(M,\mathbb{Z})$ for two compact hyperbolic 3-manifolds:"
    "\n"
    r"the Meyerhoff manifold \mfld{m003}~\cite{Meyerhoff87} and \mfld{m006},"
    "\n"
    r"both with $H_1=\mathbb{Z}/5$."
    "\n"
    r"The abelianization partitions closed geodesics into 5 homology classes;"
    "\n"
    r"we compute the spectral floor $\phi_{\rm floor}(k)=\min_{\gamma\in C_k}"
    "\n"
    r"|\phi_{\rm fold}(\gamma)|$ for each class $k\in\mathbb{Z}/5$."
    "\n\n"
    r"Our main finding is a striking asymmetry on \mfld{m006}: the five spectral"
    "\n"
    r"floors span a factor of $637$ at word lengths $\leq 7$, with class~4"
    "\n"
    r"containing a twist-suppressed loxodromic geodesic \word{AAABAB}"
    "\n"
    r"($\phi_{\rm fold}=0.0053\degr$, $|\mathrm{Tr}|=3.998$, real length $\ell=2.633$)"
    "\n"
    r"while all other classes have floors $\geq 1.61\degr$."
    "\n"
    r"This is in contrast to the Meyerhoff manifold \mfld{m003}, where the"
    "\n"
    r"five floors lie within a factor of $2$."
    "\n\n"
    r"We make the following contributions:"
    "\n"
    r"\begin{enumerate}"
    "\n"
    r"\item A systematic computation of homology-class spectral floors for"
    "\n"
    r"  \mfld{m003} and \mfld{m006} up to word length~7 (Table~\ref{tab:floors})."
    "\n"
    r"\item Identification of a twist-suppressed loxodromic element \word{AAABAB}"
    "\n"
    r"  in class~4 of \mfld{m006} with full geometric invariants"
    "\n"
    r"  ($\ell$, $\phi$, $\mathrm{Tr}$, $|\lambda|$)."
    "\n"
    r"\item A stability observation: the class~3 floor $1.611\degr$ is unchanged"
    "\n"
    r"  from word length~6 to~7, suggesting the gap boundary is stable."
    "\n"
    r"\item A falsifiable prediction: no loxodromic element with"
    "\n"
    r"  $0.005\degr < \phi_{\rm fold} < 1.61\degr$ exists in $\pi_1(\mfld{m006})$"
    "\n"
    r"  at any word length, which is testable by extending the census."
    "\n"
    r"\item A connection to the arithmetic structure of \mfld{m003}~\cite{Chinburg98}"
    "\n"
    r"  and to Standard Model flavor parameters~\cite{GentryTwist}."
    "\n"
    r"\end{enumerate}"
)

tex = tex.replace(old_intro, new_intro, 1)

# ── Fix 2: Add stability observation to results ────────────────────
old_gap = (
    r"We note also that $5\times 1.611\degr = 8.055\degr$, which is within 5.7\% of"
    "\n"
    r"the PMNS reactor angle $\theta_{13}^\nu=8.54\degr$~\cite{PDG2024}."
    "\n"
    r"This may reflect a deeper $\mathbb{Z}/5$ organisation of the spectral floors,"
    "\n"
    r"but we offer it as an observation rather than a claim."
)

new_gap = (
    r"\textbf{Stability.}"
    "\n"
    r"The class~3 floor is set by \word{ABaBAb} (length~6, $\phi_{\rm fold}=1.611\degr$)."
    "\n"
    r"Extending the census from length~6 to~7, the class~3 floor does not decrease:"
    "\n"
    r"no length-7 word in class~3 achieves $\phi_{\rm fold} < 1.611\degr$."
    "\n"
    r"This provides evidence that the gap boundary is not a finite-length artifact,"
    "\n"
    r"though we cannot exclude that longer words in class~3 could close the gap."
    "\n"
    r"Extending the census to length~10 is the most direct test."
    "\n\n"
    r"We note also that $5\times 1.611\degr = 8.055\degr$, within 5.7\% of"
    "\n"
    r"the PMNS reactor angle $\theta_{13}^\nu=8.54\degr$~\cite{PDG2024};"
    "\n"
    r"we offer this as an observation rather than a claim."
)

tex = tex.replace(old_gap, new_gap, 1)

# ── Fix 3: Remove "expected to be arithmetic" speculation ──────────
old_arith = (
    r"\mfld{m006} is expected to be arithmetic by similar criteria (small volume, torsion homology), "
    r"but its explicit invariant trace field and quaternion algebra have not been computed in the literature "
    r"to our knowledge; this is deferred to future work using the Snap package~\cite{Snap}."
)
new_arith = (
    r"The arithmeticity of \mfld{m006} has not been established in the literature to our knowledge;"
    r" its trace $\mathrm{Tr}(\word{AAABAB})=-3.998+0.0003i$ suggests the holonomy takes"
    r" values in a number field, consistent with arithmeticity, but a rigorous determination"
    r" requires computation of the invariant trace field via the Snap package~\cite{Snap},"
    r" which is deferred to future work."
)
tex = tex.replace(old_arith, new_arith, 1)

# ── Fix 4: Check for unresolved references ─────────────────────────
import re
unresolved = re.findall(r'\[\?\]|Table \?\?|Fig\.\s*\?\?', tex)
if unresolved:
    print(f"WARNING: Found unresolved references: {unresolved}")
else:
    print("No unresolved [?] references found.")

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("All fixes applied.")
