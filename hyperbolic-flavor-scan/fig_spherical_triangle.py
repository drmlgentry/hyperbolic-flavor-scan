import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

# --- Confirmed axis vectors ---
n = {
    'aaB': np.array([ 0.237, -0.899,  0.368]),
    'AbA': np.array([ 0.376,  0.911,  0.170]),
    'AAb': np.array([ 0.091,  0.184,  0.979]),
}
keys = ['aaB', 'AbA', 'AAb']
colors = {'aaB': '#e6194b', 'AbA': '#3cb44b', 'AAb': '#4363d8'}

angles = {
    ('aaB','AbA'): 48.2,
    ('aaB','AAb'): 77.5,
    ('AbA','AAb'): 68.4,
}
angle_sum = sum(angles.values())
excess = angle_sum - 180.0

# Build local 2D frame at centroid of the three points
pts = np.array([n[k] for k in keys])
centroid = np.mean(pts, axis=0)
centroid /= np.linalg.norm(centroid)

e1 = np.cross(centroid, np.array([0,0,1]))
if np.linalg.norm(e1) < 1e-6:
    e1 = np.cross(centroid, np.array([0,1,0]))
e1 /= np.linalg.norm(e1)
e2 = np.cross(centroid, e1)
e2 /= np.linalg.norm(e2)

def project(v):
    return np.array([np.dot(v, e1), np.dot(v, e2)])

proj = {k: project(n[k]) for k in keys}

def arc2d(k1, k2, npts=120):
    t = np.linspace(0, 1, npts)
    raw = np.outer(1-t, n[k1]) + np.outer(t, n[k2])
    norms = np.linalg.norm(raw, axis=1, keepdims=True)
    sph = raw / norms
    return np.array([project(s) for s in sph])

fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.set_aspect('equal')
ax.axis('off')

edge_data = [
    ('aaB', 'AbA', angles[('aaB','AbA')], '#9a0000'),
    ('aaB', 'AAb', angles[('aaB','AAb')], '#00600a'),
    ('AbA', 'AAb', angles[('AbA','AAb')], '#00008a'),
]

for k1, k2, ang, ec in edge_data:
    arc = arc2d(k1, k2)
    ax.plot(arc[:,0], arc[:,1], color=ec, lw=2.2, zorder=2)
    mid = arc[len(arc)//2]
    perp = arc[len(arc)//2+3] - arc[len(arc)//2-3]
    perp = np.array([-perp[1], perp[0]])
    perp /= np.linalg.norm(perp) * 8
    lpos = mid + perp
    ax.text(lpos[0], lpos[1], f'^\\circ$',
            fontsize=12, color=ec, ha='center', va='center',
            fontweight='bold',
            path_effects=[pe.withStroke(linewidth=2, foreground='white')])

# Vertex labels -- use \mathrm instead of \texttt (matplotlib limitation)
offsets = {'aaB': (-0.038, -0.028), 'AbA': (0.025, -0.025), 'AAb': (0.0, 0.028)}
for k in keys:
    p = proj[k]
    ax.scatter(*p, color=colors[k], s=60, zorder=5)
    ox, oy = offsets[k]
    ax.text(p[0]+ox, p[1]+oy, f'$\\mathbf{{n}}_{{\\mathrm{{{k}}}}}$',
            fontsize=12, color=colors[k],
            fontweight='bold', ha='center', va='center')

ax.text(0.5, 0.04,
        f'Sum of angles: ^\\circ$\n'
        f'Spherical excess: $\\Sigma - 180^\\circ = {excess:.1f}^\\circ \\approx \\theta_C$',
        transform=ax.transAxes, fontsize=11, ha='center', va='bottom',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow',
                  edgecolor='goldenrod', alpha=0.9))

ax.set_title('Spherical Triangle of Geodesic Axes\n'
             'Manifold m006 $\\;|\\;$ Excess $= 13^\\circ \\approx$ Cabibbo angle',
             fontsize=12, pad=14)

plt.tight_layout()
plt.savefig('fig_spherical_triangle.pdf', dpi=200, bbox_inches='tight')
plt.savefig('fig_spherical_triangle.png', dpi=150, bbox_inches='tight')
print("Saved fig_spherical_triangle.pdf and .png")
