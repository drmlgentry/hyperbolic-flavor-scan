import snappy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

phi = (1 + np.sqrt(5)) / 2
log_phi = np.log(phi)
gold  = '#e8c96a'
teal  = '#5ec8c8'
white = '#f0f0f0'
gold2 = '#c8a84b'
bg    = '#0a0e1a'
axbg  = '#0e1528'

def get_lengths(M, cutoff=5.0):
    return np.array(sorted([float(x.length.real()) for x in M.length_spectrum(cutoff)]))

def build_kernel(L, s):
    d = L[:,None] - L[None,:]
    return np.exp(-d*d/(2*s*s))

def eff_rank(K):
    ev = np.linalg.eigvalsh(K)
    ev = np.maximum(ev, 1e-12)
    p = ev/ev.sum()
    return np.exp(-np.sum(p*np.log(p)))

def energy_frac(K, n):
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    return ev[:n].sum()/ev.sum()

print("Loading manifolds...", flush=True)
M1 = snappy.Manifold('m003'); M1.dehn_fill((-2,3))
M2 = snappy.Manifold('m006'); M2.dehn_fill((-5,2))
L1 = get_lengths(M1, 5.0)
L2 = get_lengths(M2, 5.0)
print("PMNS: %d  CKM: %d geodesics" % (len(L1),len(L2)), flush=True)

sigmas = np.linspace(0.20, 1.30, 80)
print("Computing...", flush=True)
R1 = np.array([eff_rank(build_kernel(L1,s)) for s in sigmas])
R2 = np.array([eff_rank(build_kernel(L2,s)) for s in sigmas])
E3_1 = np.array([energy_frac(build_kernel(L1,s),3) for s in sigmas])
E3_2 = np.array([energy_frac(build_kernel(L2,s),3) for s in sigmas])
E2_1 = np.array([energy_frac(build_kernel(L1,s),2) for s in sigmas])
E2_2 = np.array([energy_frac(build_kernel(L2,s),2) for s in sigmas])
print("Plotting...", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.patch.set_facecolor(bg)
for ax in axes:
    ax.set_facecolor(axbg)
    ax.tick_params(colors=gold2, labelsize=11)
    for sp in ax.spines.values():
        sp.set_edgecolor(gold2); sp.set_linewidth(0.8)

win_lo, win_hi = np.sqrt(1/6), np.sqrt(1/3)

ax = axes[0]
ax.plot(sigmas, R1, color=gold, linewidth=2.4, label='PMNS = m003(-2,3)')
ax.plot(sigmas, R2, color=teal, linewidth=2.4, linestyle='--', label='CKM  = m006(-5,2)')
ax.axvspan(win_lo, win_hi, alpha=0.13, color=teal, label='3-gen window')
ax.axvline(log_phi, color=gold, linestyle='-', linewidth=1.3, alpha=0.85)
ax.axvline(1.5*log_phi, color=teal, linestyle='--', linewidth=1.3, alpha=0.85)
ax.axhline(3, color=white, linestyle=':', linewidth=0.8, alpha=0.5)
R1_lp = float(R1[np.argmin(np.abs(sigmas-log_phi))])
R2_lp = float(R2[np.argmin(np.abs(sigmas-log_phi))])
ax.annotate('$R_{\\rm eff}=%.2f$ (PMNS)' % R1_lp,
            xy=(log_phi, R1_lp), xytext=(log_phi+0.12, R1_lp+0.3),
            color=gold, fontsize=10,
            arrowprops=dict(arrowstyle='->', color=gold, lw=1.2))
ax.annotate('$R_{\\rm eff}=%.2f$ (CKM)' % R2_lp,
            xy=(log_phi, R2_lp), xytext=(log_phi+0.12, R2_lp-0.6),
            color=teal, fontsize=10,
            arrowprops=dict(arrowstyle='->', color=teal, lw=1.2))
ax.text(log_phi-0.03, 1.6, '$\\log\\phi$', color=gold, fontsize=9, ha='right', rotation=90)
ax.text(1.5*log_phi+0.02, 1.6, '$(3/2)\\log\\phi$', color=teal, fontsize=9, rotation=90)
ax.text(win_lo+0.01, 5.6, '$\\sqrt{1/6}$', color=white, fontsize=8)
ax.text(win_hi+0.01, 5.6, '$\\sqrt{1/3}$', color=white, fontsize=8)
ax.text(1.22, 3.08, 'R=3', color=white, fontsize=8)
ax.set_xlabel('$\\sigma$', color=white, fontsize=13)
ax.set_ylabel('Effective rank $R_{\\rm eff}$', color=white, fontsize=12)
ax.set_title('Geodesic length kernel: soft dimensional transition', color=white, fontsize=11, pad=8)
ax.set_xlim(0.20, 1.30); ax.set_ylim(1.5, 6.5)
ax.legend(fontsize=9, facecolor=axbg, edgecolor=gold2, labelcolor=white, loc='upper right')
ax.grid(alpha=0.15, color=gold2)

ax = axes[1]
ax.plot(sigmas, E3_1*100, color=gold, linewidth=2.4, label='Top 3 modes (PMNS)')
ax.plot(sigmas, E3_2*100, color=teal, linewidth=2.4, linestyle='--', label='Top 3 modes (CKM)')
ax.plot(sigmas, E2_1*100, color=gold, linewidth=1.2, linestyle=':', alpha=0.6, label='Top 2 modes (PMNS)')
ax.plot(sigmas, E2_2*100, color=teal, linewidth=1.2, linestyle=':', alpha=0.6, label='Top 2 modes (CKM)')
ax.axvspan(win_lo, win_hi, alpha=0.13, color=teal)
ax.axvline(log_phi, color=gold, linestyle='-', linewidth=1.3, alpha=0.85)
ax.axvline(1.5*log_phi, color=teal, linestyle='--', linewidth=1.3, alpha=0.85)
E3_1_lp = float(E3_1[np.argmin(np.abs(sigmas-log_phi))])*100
E3_2_lp = float(E3_2[np.argmin(np.abs(sigmas-log_phi))])*100
ax.axhline(E3_1_lp, color=gold, linestyle=':', linewidth=0.8, alpha=0.5)
ax.axhline(E3_2_lp, color=teal, linestyle=':', linewidth=0.8, alpha=0.5)
ax.annotate('%.1f%% PMNS' % E3_1_lp,
            xy=(log_phi, E3_1_lp), xytext=(log_phi+0.10, E3_1_lp-4),
            color=gold, fontsize=10,
            arrowprops=dict(arrowstyle='->', color=gold, lw=1.2))
ax.annotate('%.1f%% CKM' % E3_2_lp,
            xy=(log_phi, E3_2_lp), xytext=(log_phi+0.10, E3_2_lp+2),
            color=teal, fontsize=10,
            arrowprops=dict(arrowstyle='->', color=teal, lw=1.2))
ax.text(log_phi-0.03, 62, '$\\log\\phi$', color=gold, fontsize=9, ha='right', rotation=90)
ax.text(1.5*log_phi+0.02, 62, '$(3/2)\\log\\phi$', color=teal, fontsize=9, rotation=90)
ax.set_xlabel('$\\sigma$', color=white, fontsize=13)
ax.set_ylabel('Spectral energy fraction (%)', color=white, fontsize=12)
ax.set_title('Energy in top modes vs bandwidth', color=white, fontsize=11, pad=8)
ax.set_xlim(0.20, 1.30); ax.set_ylim(60, 100)
ax.legend(fontsize=8, facecolor=axbg, edgecolor=gold2, labelcolor=white, loc='lower right')
ax.grid(alpha=0.15, color=gold2)

fig.suptitle('Soft dimensional transition at $\\sigma=\\log\\phi$: geodesic spectra of $M_{\\rm PMNS}$ and $M_{\\rm CKM}$',
             color=white, fontsize=12, y=1.01)
plt.subplots_adjust(wspace=0.32)

plt.savefig('fig_kernel_rank.png', dpi=200, bbox_inches='tight', facecolor=bg)
plt.savefig('fig_kernel_rank.pdf', dpi=200, bbox_inches='tight', facecolor=bg)
plt.close()
print("Saved fig_kernel_rank.png and fig_kernel_rank.pdf", flush=True)
