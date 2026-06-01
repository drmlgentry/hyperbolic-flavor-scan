path = r"C:\dev\framework\papers\flavor-mixing\gentry-flavor-mixing.tex"
with open(path, encoding="utf-8") as f:
    tex = f.read()

old_remark = (
    r"\begin{remark}[Poisson kernel realization]" + "\n"
    r"An explicit family satisfying Definition~\ref{def:localized} can be" + "\n"
    r"constructed using the Poisson kernel for $\Hthree$." + "\n"
    r"For $\xi\in\Sphere$ and a basepoint $o\in\Hthree$, set" + "\n"
    r"$\ket{\xi}$ to be the $L^2(\Sphere)$ function proportional to the" + "\n"
    r"Poisson kernel $P(\cdot,\xi)$, normalized in $L^2(\Sphere)$." + "\n"
    r"The resulting family is continuous in $\xi$, and" + "\n"
    r"$\braket{\xi|\zeta}=1$ if and only if $\xi=\zeta$~\cite{ratcliffe}." + "\n"
    r"\end{remark}"
)
new_remark = (
    r"\begin{remark}[Poisson kernel realization]" + "\n"
    r"An explicit family satisfying Definition~\ref{def:localized} can be" + "\n"
    r"constructed using the Poisson kernel for $\Hthree$." + "\n"
    r"For $\xi\in\Sphere$ and a basepoint $o\in\Hthree$, set" + "\n"
    r"$\ket{\xi}$ to be the $L^2(\Sphere)$ function proportional to the" + "\n"
    r"Poisson kernel $P(\cdot,\xi)$, normalized in $L^2(\Sphere)$." + "\n"
    r"The resulting family is continuous in $\xi$, satisfies" + "\n"
    r"$\braket{\xi|\xi}=1$, and $|\braket{\xi|\zeta}|<1$ for $\xi\neq\zeta$~\cite{ratcliffe}." + "\n"
    r"In particular, distinct boundary directions produce linearly independent" + "\n"
    r"vectors, which is the property required by Definition~\ref{def:localized}." + "\n"
    r"\end{remark}"
)
tex = tex.replace(old_remark, new_remark, 1)

with open(path, "w", encoding="utf-8") as f:
    f.write(tex)
print("Done." if new_remark[:30] in open(path, encoding="utf-8").read() else "FAILED.")
