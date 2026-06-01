import re, sys

def fix_paper(src_path, out_path, title, extra_theorems=True):
    with open(src_path, encoding="utf-8-sig") as f:
        raw = f.read()

    # Find \begin{document} and split there
    doc_idx = raw.find(r"\begin{document}")
    preamble = raw[:doc_idx]
    body = raw[doc_idx:]

    # Build clean preamble from scratch — keep only \usepackage and \newcommand lines
    clean_pre = [r"\documentclass[11pt,a4paper]{article}", ""]
    for line in preamble.splitlines():
        s = line.strip()
        if s.startswith(r"\usepackage") and not any(x in s for x in ["svglov3","spbasic","svjour3"]):
            clean_pre.append(line)
        elif s.startswith(r"\newcommand") or s.startswith(r"\renewcommand"):
            clean_pre.append(line)
        elif s.startswith(r"\newcolumntype"):
            clean_pre.append(line)

    clean_pre += [
        "",
        r"\usepackage{authblk}",
        r"\usepackage[margin=2.5cm]{geometry}",
        r"\usepackage{amsthm}",
        r"\newtheorem{proposition}{Proposition}",
        r"\newtheorem{theorem}{Theorem}",
        r"\newtheorem{lemma}{Lemma}",
        r"\newtheorem{corollary}{Corollary}",
        "",
        r"\title{" + title + "}",
        r"\author[1]{Marvin L. Gentry}",
        r"\affil[1]{Independent Researcher, Seattle, WA, USA. \texttt{drmlgentry@protonmail.com}. ORCID: 0009-0006-4550-2663}",
        r"\date{\today}",
        "",
    ]

    # Clean body: remove svjour3 author commands
    skip = [r"\author{", r"\author[", r"\institute{", r"\email{",
            r"\titlerunning{", r"\authorrunning{", r"\title{", r"\maketitle"]
    clean_body = []
    for line in body.splitlines(keepends=True):
        if any(line.strip().startswith(p) for p in skip):
            continue
        clean_body.append(line)
    body = "".join(clean_body)

    body = body.replace(r"\begin{document}", r"\begin{document}" + "\n" + r"\maketitle", 1)
    body = body.replace(r"\bibliographystyle{spbasic}", r"\bibliographystyle{unsrt}")
    body = body.replace(r"\begin{acknowledgments}", r"\section*{Acknowledgements}")
    body = body.replace(r"\end{acknowledgments}", "")
    body = body.replace(r"\begin{acknowledgement}", r"\section*{Acknowledgements}")
    body = body.replace(r"\end{acknowledgement}", "")

    result = "\n".join(clean_pre) + "\n" + body
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(result)
    print(f"Written: {out_path} ({len(result.splitlines())} lines)")

fix_paper(
    r"C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-hyperbolic-flavor-cp-epjc.tex",
    r"C:\dev\framework\papers\hyperbolic-flavor-cp\gentry-cp-rip.tex",
    r"CP Violation from A-Factor Twist Angles of Hyperbolic Holonomy"
)

fix_paper(
    r"C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-hyperbolic-flavor-twist-epjc.tex",
    r"C:\dev\framework\papers\hyperbolic-flavor-twist\gentry-twist-rip.tex",
    r"Twist Angle Spectrum of Hyperbolic Holonomy: Encoding of Standard Model Flavor Parameters"
)
