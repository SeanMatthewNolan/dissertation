import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

alg_set = gp.load_sol_set('downrange_alg.data')
dif_set = gp.load_sol_set('downrange_dif.data')
# com_sol = gp.load_sol_set('downrange_com.data')[-1]

com_sol = dif_set[-1]

gam_0s_alg = np.rad2deg(np.array([sol.x[3, 0] for sol in alg_set]))
gam_0s_dif = np.rad2deg(np.array([sol.x[3, 0] for sol in dif_set]))

cmap = plt.colormaps['viridis']
norm = Normalize(min(gam_0s_alg), max(gam_0s_alg))

fig = plt.figure(figsize=(6.5, 3))
title = fig.suptitle(r'$H_{\bm{uu}}$ at High Altitudes')

gs = gridspec.GridSpec(2, 3, height_ratios=[1, 0.1])

ax_alt = fig.add_subplot(gs[0, 0])
ax_h_uu = fig.add_subplot(gs[0, 1])
ax_h_uu_alt = fig.add_subplot(gs[0, 2])

ax_cb = fig.add_subplot(gs[1, :])

for sol in alg_set:
    color = cmap(norm(np.rad2deg(sol.x[3, 0])))

    ax_alt.plot(sol.t * 180 / np.pi, sol.x[0, :] / 1e5, color=color)
    ax_h_uu.plot(sol.t, sol.aux['H_uu'][:, 0, 0], color=color)
    ax_h_uu_alt.plot(sol.x[0, :] / 1e5, sol.aux['H_uu'][:, 0, 0], color=color)

ax_alt.plot(com_sol.t * 180 / np.pi, com_sol.x[0, :] / 1e5, label='Last Diff. Sol.', color='C3')
ax_h_uu.plot(com_sol.t, com_sol.aux['H_uu'][:, 0, 0], label='Last Diff. Sol.', color='C3')
ax_h_uu_alt.plot(com_sol.x[0, :] / 1e5, com_sol.aux['H_uu'][:, 0, 0], label='Last Diff. Sol.', color='C3')

fig.colorbar(plt.cm.ScalarMappable(norm=norm), cmap=cmap,
             cax=ax_cb, label=r'Initial FPA $\gamma_0$, [deg]',
             orientation='horizontal')

# ax_alt.set_title('Algebraic Control Law: Path')
ax_alt.set_xlabel(r'Time, $t$ [s]')
ax_alt.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')
ax_alt.legend(frameon=False, loc='upper center')
ax_alt.set_ylim(0, 25)

# ax_h_uu.set_title('Algebraic Control Law: AOA')
ax_h_uu.set_xlabel(r'Time $t$ [s]')
ax_h_uu.set_ylabel(r'$H_{\bm{uu}}$ [1/rad/s]')
ax_h_uu.set_yscale('log')

# ax_h_uu_alt.set_title('Differential Control Law: Path')
ax_h_uu_alt.set_xlabel(r'Altitude, $h$ [$\times10^5$ ft]')
ax_h_uu_alt.set_ylabel(r'$H_{\bm{uu}}$ [1/rad/s]')
ax_h_uu_alt.set_yscale('log')

fig.subplots_adjust(
        top=0.851,
        bottom=0.188,
        left=0.091,
        right=0.977,
        hspace=0.792,
        wspace=0.562
)

min_h_uu = min(com_sol.aux['H_uu'][:, 0, 0])
print(max(gam_0s_dif))
print(min_h_uu)
print(f'{1/min_h_uu:e}')

plt.show()
