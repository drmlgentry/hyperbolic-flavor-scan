import subprocess, os
# Try converting PNG to PDF via matplotlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

img = mpimg.imread(r"C:\dev\hyperbolic-flavor-scan\figures\decay_multiclass_m006_clean.png")
fig, ax = plt.subplots(figsize=(8, 6))
ax.imshow(img)
ax.axis("off")
fig.savefig(r"C:\dev\framework\papers\hyperbolic-flavor-torsion\fig_floor_decay.pdf",
            format="pdf", dpi=150, bbox_inches="tight", pad_inches=0)
plt.close()
print("Converted to PDF.")
