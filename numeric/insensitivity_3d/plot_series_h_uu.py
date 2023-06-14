import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

alg_set = gp.load_sol_set('crossrange.data')

gam_0s_alg = np.rad2deg(np.array([sol.x[4, 0] for sol in alg_set]))

cmap = plt.colormaps['viridis']
norm = Normalize(min(gam_0s_alg), max(gam_0s_alg))

fig = plt.figure(figsize=(6.5, 3))
title = fig.suptitle(r'Eigenvalues of $H_{\bm{uu}}$ at High Altitudes:' + '\nSpace Shuttle Crossrange Problem')

gs = gridspec.GridSpec(2, 3, height_ratios=[1, 0.1])

ax_alt = fig.add_subplot(gs[0, 0])
ax_h_uu = fig.add_subplot(gs[0, 1])
ax_h_uu_alt = fig.add_subplot(gs[0, 2])

ax_cb = fig.add_subplot(gs[1, :])

min_eig_1s = []
min_eig_2s = []

for sol in alg_set:
    color = cmap(norm(np.rad2deg(sol.x[4, 0])))

    eig_h_uu = np.array([np.linalg.eig(h_uu)[0] for h_uu in sol.aux['H_uu']])

    min_eig_1s.append(min(eig_h_uu[:, 0]))
    min_eig_2s.append(min(eig_h_uu[:, 1]))

    ax_alt.plot(sol.t * 180 / np.pi, sol.x[0, :] / 1e5, color=color)
    ax_h_uu.plot(sol.t, eig_h_uu[:, 0], color=color)
    ax_h_uu_alt.plot(sol.t, eig_h_uu[:, 1], color=color)

fig.colorbar(plt.cm.ScalarMappable(norm=norm), cmap=cmap,
             cax=ax_cb, label=r'Initial FPA $\gamma_0$, [deg]',
             orientation='horizontal')

# ax_alt.set_title('Algebraic Control Law: Path')
ax_alt.set_xlabel(r'Time, $t$ [s]')
ax_alt.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

# ax_h_uu.set_title('Algebraic Control Law: AOA')
ax_h_uu.set_xlabel(r'Time $t$ [s]')
ax_h_uu.set_ylabel(r'$\text{eig}_{1}\left({H_{\bm{uu}}}\right)$')
ax_h_uu.set_yscale('log')

# ax_h_uu_alt.set_title('Differential Control Law: Path')
ax_h_uu_alt.set_xlabel(r'Time $t$ [s]')
ax_h_uu_alt.set_ylabel(r'$\text{eig}_{2}\left({H_{\bm{uu}}}\right)$')
ax_h_uu_alt.set_yscale('log')

fig.subplots_adjust(
        top=0.794,
        bottom=0.188,
        left=0.08,
        right=0.977,
        hspace=0.9,
        wspace=0.572
)

print(gam_0s_alg[-1])
print(min(min_eig_1s))
print(min(min_eig_2s))

plt.show()
