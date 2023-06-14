import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

# sol_set = gp.load_sol_set('downrange_alg_10.data')
sol_set = gp.load_sol_set('downrange_dif_mit.data')
# com_sol = gp.load_sol_set('downrange_com.data')[-1]

gam_0s = np.rad2deg(np.array([sol.x[3, 0] for sol in sol_set]))

cmap = plt.colormaps['viridis']
norm = Normalize(min(gam_0s), round(max(gam_0s)))

fig = plt.figure(figsize=(3.5, 2.9))
title = fig.suptitle(r'Mitigation of Control Insensitivity')

gs = gridspec.GridSpec(2, 2, width_ratios=[1, 0.05])

ax_alpha = fig.add_subplot(gs[0, 0])
ax_h_uu = fig.add_subplot(gs[1, 0])

ax_cb = fig.add_subplot(gs[:, 1])

for sol in sol_set:
    color = cmap(norm(np.rad2deg(sol.x[3, 0])))

    ax_alpha.plot(sol.t, np.rad2deg(sol.u[0, :]), color=color)
    ax_h_uu.plot(sol.t, sol.aux['H_uu'][:, 0, 0], color=color)

fig.colorbar(plt.cm.ScalarMappable(norm=norm), cmap=cmap,
             cax=ax_cb, label=r'Initial FPA $\gamma_0$, [deg]',
             orientation='vertical')


ax_alpha.set_xlabel(r'Time, $t$ [s]')
ax_alpha.set_ylabel(r'AOA, $\alpha$ [deg]')

# ax_h_uu.set_title('Algebraic Control Law: AOA')
ax_h_uu.set_xlabel(r'Time $t$ [s]')
ax_h_uu.set_ylabel(r'$H_{\bm{uu}}$ [1/rad/s]')
ax_h_uu.set_yscale('log')

fig.subplots_adjust(
        top=0.829,
        bottom=0.195,
        left=0.206,
        right=0.834,
        hspace=0.884,
        wspace=0.146
)

min_h_uu = min(sol_set[-1].aux['H_uu'][:, 0, 0])
print(max(gam_0s))
print(min_h_uu)
print(f'{1/min_h_uu:e}')

plt.show()
