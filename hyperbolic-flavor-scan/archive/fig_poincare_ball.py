import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# --- Confirmed axis vectors (Pauli/Doc-8 method) ---
axes = {
    'aaB': np.array([ 0.237, -0.899,  0.368]),
    'AbA': np.array([ 0.376,  0.911,  0.170]),
    'AAb': np.array([ 0.091,  0.184,  0.979]),
}
colors = {'aaB': '#e6194b', 'AbA': '#3cb44b', 'AAb': '#4363d8'}
labels = {'aaB': r'$\gamma_1$ (aaB)', 'AbA': r'$\gamma_2$ (AbA)', 'AAb': r'$\gamma_3$ (AAb)'}

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')

# Draw transparent unit sphere
u = np.linspace(0, 2*np.pi, 60)
v = np.linspace(0, np.pi, 40)
xs = np.outer(np.cos(u), np.sin(v))
ys = np.outer(np.sin(u), np.sin(v))
zs = np.outer(np.ones_like(u), np.cos(v))
ax.plot_surface(xs, ys, zs, alpha=0.07, color='steelblue', linewidth=0)

# Draw faint latitude/longitude lines
for phi in np.linspace(0, np.pi, 7):
    x = np.sin(phi)*np.cos(u)
    y = np.sin(phi)*np.sin(u)
    z = np.cos(phi)*np.ones_like(u)
    ax.plot(x, y, z, color='gray', lw=0.3, alpha=0.4)
for th in np.linspace(0, 2*np.pi, 13):
    vv = np.linspace(0, np.pi, 40)
    x = np.sin(vv)*np.cos(th)
    y = np.sin(vv)*np.sin(th)
    z = np.cos(vv)
    ax.plot(x, y, z, color='gray', lw=0.3, alpha=0.4)

# Draw geodesic arcs between each pair
def geodesic_arc(n1, n2, npts=80):
    t = np.linspace(0, 1, npts)
    pts = np.outer(1-t, n1) + np.outer(t, n2)
    norms = np.linalg.norm(pts, axis=1, keepdims=True)
    return pts / norms

pairs = [('aaB','AbA', r'.2^\circ$'),
         ('aaB','AAb', r'.5^\circ$'),
         ('AbA','AAb', r'.4^\circ$')]

arc_colors = ['#9a0000', '#00600a', '#00008a']
for (k1, k2, ang), ac in zip(pairs, arc_colors):
    arc = geodesic_arc(axes[k1], axes[k2])
    ax.plot(arc[:,0], arc[:,1], arc[:,2], color=ac, lw=1.8, alpha=0.85)
    mid = arc[len(arc)//2]
    mid_n = mid / np.linalg.norm(mid) * 1.13
    ax.text(mid_n[0], mid_n[1], mid_n[2], ang, fontsize=9,
            ha='center', va='center', color=ac)

# Draw axis arrows and labels
for name, n in axes.items():
    c = colors[name]
    ax.quiver(0, 0, 0, n[0], n[1], n[2],
              length=1.0, color=c, linewidth=2.2,
              arrow_length_ratio=0.12)
    tip = n * 1.18
    ax.text(tip[0], tip[1], tip[2], labels[name],
            fontsize=10, color=c, fontweight='bold', ha='center')

# Coordinate axes (thin)
for d, lbl in zip(np.eye(3), ['x','y','z']):
    ax.quiver(0,0,0,*d*1.25, color='black', lw=0.8, arrow_length_ratio=0.08, alpha=0.4)
    ax.text(*(d*1.32), lbl, fontsize=8, color='gray', ha='center')

ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_zlim(-1.2, 1.2)
ax.set_box_aspect([1,1,1])
ax.set_axis_off()
ax.view_init(elev=22, azim=38)

ax.set_title('Geodesic Axes on $\\partial\\mathbb{H}^3 \\cong S^2$\n'
             'Manifold m006, words aaB / AbA / AAb',
             fontsize=11, pad=10)

plt.tight_layout()
plt.savefig('fig_poincare_ball.pdf', dpi=200, bbox_inches='tight')
plt.savefig('fig_poincare_ball.png', dpi=150, bbox_inches='tight')
print("Saved fig_poincare_ball.pdf and .png")
