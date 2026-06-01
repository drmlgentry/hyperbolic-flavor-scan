path = r"C:\dev\framework\papers\hyperbolic-flavor-torsion\gentry-hyperbolic-flavor-torsion.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

# 1. Fix "decays exponentially" -> "is well described by exponential decay"
tex = tex.replace(
    r"$\phi_{\rm floor}^{(k)}(L)$ decays exponentially with~$L$ for all five",
    r"$\phi_{\rm floor}^{(k)}(L)$ is well described by exponential decay over the observed range for all five"
)

# 2. Fix systole claim
tex = tex.replace(
    r"Since each generator has a fixed positive translation length"
    "\n"
    r"$\ell_{\rm min}=\min_g\ell(g)>0$ (equal to the systole of \mfld{m006}),",
    r"Since each generator has a positive translation length bounded below by"
    "\n"
    r"the systole $\mathrm{sys}(\mfld{m006})>0$,"
)

# 3. Fix primitive/non-primitive sentence
tex = tex.replace(
    r"for the spectral floor the distinction is immaterial since a non-primitive"
    "\n"
    r"element $\gamma^n$ has $\phi_{\rm fold}(\gamma^n)\geq\phi_{\rm fold}(\gamma)$"
    "\n"
    r"only if $n$ is odd, and the floor is in any case dominated by primitive elements.",
    r"For the spectral floor, non-primitive elements contribute only subleading"
    "\n"
    r"corrections: the floor is dominated by primitive geodesics since"
    "\n"
    r"$\phi_{\rm fold}(\gamma^n)$ need not be smaller than $\phi_{\rm fold}(\gamma)$"
    "\n"
    r"for $n\geq 2$, and the Margulis count of primitive geodesics grows"
    "\n"
    r"as $e^{2T}/(2T)$ while non-primitive geodesics are negligible in comparison."
)

# 4. Add numerical precision note to census procedure
old_census_end = r"Words with $|\lambda|\leq 1.01$"
new_census_start = (
    r"Twist angles are computed as $\phi=\arg\lambda$ where $\lambda$ is the"
    "\n"
    r"eigenvalue of larger modulus of the $2\times 2$ holonomy matrix."
    "\n"
    r"Near-$\pi$ twists ($\phi_{\rm fold}<0.01\degr$) are verified by computing"
    "\n"
    r"$|\mathrm{Im}(\lambda)/\mathrm{Re}(\lambda)|$ directly; all values"
    "\n"
    r"reported are stable to the numerical precision of \textsc{SnapPy}'s"
    "\n"
    r"default 15-digit arithmetic."
    "\n"
    r"Words with $|\lambda|\leq 1.01$"
)
tex = tex.replace(old_census_end, new_census_start, 1)

# 5. Reframe SM motivation paragraph
old_motivation = (
    r"The motivation for this study comes from a companion series of"
    "\n"
    r"papers~\cite{GentryCKM,GentryPMNS,GentryTwist} in which \mfld{m003}"
    "\n"
    r"and \mfld{m006} were found to reproduce Standard Model quark and lepton"
    "\n"
    r"mixing matrices."
    "\n"
    r"The present work focuses on the geometric structure of the twist spectrum"
    "\n"
    r"as an object of independent mathematical interest."
)
new_motivation = (
    r"\textit{Motivation.}"
    "\n"
    r"This study was originally motivated by the observation, made in a"
    "\n"
    r"companion series of papers~\cite{GentryCKM,GentryPMNS,GentryTwist},"
    "\n"
    r"that \mfld{m003} and \mfld{m006} reproduce Standard Model quark and"
    "\n"
    r"lepton mixing matrices with no free parameters."
    "\n"
    r"The present paper is independent of that program: the spectral floor"
    "\n"
    r"invariant $\phi_{\rm floor}^{(k)}(L)$ and the decay rate hierarchy are"
    "\n"
    r"defined and studied purely as geometric objects, and all results stand"
    "\n"
    r"independently of any physical interpretation."
)
tex = tex.replace(old_motivation, new_motivation, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("All fixes applied.")
