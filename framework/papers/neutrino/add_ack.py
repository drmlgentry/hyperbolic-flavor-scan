path = "gentry-neutrino.tex"
with open(path, encoding="utf-8") as f:
    src = f.read()

ack = r"""
\section*{Acknowledgements}
The author used DeepSeek and Claude (Anthropic) as AI assistants during the
preparation of this manuscript, including for \LaTeX\ editing and
computational script development. All scientific content,
mathematical results, and conclusions are the sole responsibility
of the author.

"""

src = src.replace(r"\section*{Data Availability}", ack + r"\section*{Data Availability}")

with open(path, "w", encoding="utf-8") as f:
    f.write(src)
print("Done")
