import re

for fname in ["twist_census.py", "selberg_setup.py"]:
    path = r"C:\dev\hyperbolic-flavor-scan\\" + fname
    with open(path, "r", encoding="utf-8") as f:
        txt = f.read()
    # Replace M.volume() with float(M.volume()) everywhere
    txt = txt.replace("M.volume()", "float(M.volume())")
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)
    print("Fixed: " + fname)
