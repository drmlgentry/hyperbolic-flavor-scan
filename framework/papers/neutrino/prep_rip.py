path = "gentry-neutrino.tex"
with open(path, encoding="utf-8") as f:
    src = f.read()

# Update Gentry:CORE reference to include SSRN links
src = src.replace(
    r"""  M.~L.~Gentry,
  ``Geometric Unification of Flavor: Masses, Mixing, and CP
  from the Golden Ratio Lattice,''
  preprint (2026).""",
    r"""  M.~L.~Gentry,
  ``Geometric Unification of Flavor: Masses, Mixing, and CP
  from the Golden Ratio Lattice,''
  preprint (2026). See also companion papers SSRN:6583549--6583553."""
)

"""
src = src.replace(r"\section*{Data Availability}", ack + r"\section*{Data Availability}")

with open(path, "w", encoding="utf-8") as f:
    f.write(src)
print("Done")
