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
com_sol = gp.load_sol_set('downrange_com.data')[-1]
# dif_set = gp.load_sol_set('downrange_dif_mit.data')

gam_0s_alg = np.rad2deg(np.array([sol.x[3, 0] for sol in alg_set]))
gam_0s_dif = np.rad2deg(np.array([sol.x[3, 0] for sol in dif_set]))

cmap = plt.colormaps['viridis']
norm = Normalize(min(gam_0s_alg), max(gam_0s_alg))

fig = plt.figure(figsize=(6.5, 7.25))
title = fig.suptitle('Failure of Differential Control Law at High Altitudes')

gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 0.1])

ax_path_alg = fig.add_subplot(gs[0, 0])
# ax_hv_alg = fig.add_subplot(gs[0, 1])
ax_alpha_alg = fig.add_subplot(gs[0, 1])

ax_path_dif = fig.add_subplot(gs[1, 0])
# ax_hv_dif = fig.add_subplot(gs[1, 1])
ax_alpha_dif = fig.add_subplot(gs[1, 1])

ax_path_com = fig.add_subplot(gs[2, 0])
ax_alpha_com = fig.add_subplot(gs[2, 1])

ax_cb = fig.add_subplot(gs[3, :])


for sol in alg_set:
    color = cmap(norm(np.rad2deg(sol.x[3, 0])))

    ax_path_alg.plot(sol.x[1, :] * 180 / np.pi, sol.x[0, :] / 1e5, color=color)
    # ax_hv_alg.plot(sol.x[2, :] / 1e3, sol.x[0, :] / 1e5, color=color)
    ax_alpha_alg.plot(sol.t, sol.u[0, :] * 180 / np.pi, color=color)

for sol in dif_set:
    color = cmap(norm(np.rad2deg(sol.x[3, 0])))

    ax_path_dif.plot(sol.x[1, :] * 180 / np.pi, sol.x[0, :] / 1e5, color=color)
    # ax_hv_dif.plot(sol.x[2, :] / 1e3, sol.x[0, :] / 1e5, color=color)
    ax_alpha_dif.plot(sol.t, sol.u[0, :] * 180 / np.pi, color=color)

ax_path_com.plot(com_sol.x[1, :] * 180 / np.pi, com_sol.x[0, :] / 1e5, label='Algebraic Control Law')
ax_path_com.plot(dif_set[-1].x[1, :] * 180 / np.pi, dif_set[-1].x[0, :] / 1e5, '--', label='Differential Control Law')
ax_path_com.legend()

ax_alpha_com.plot(com_sol.t, com_sol.u[0, :] * 180 / np.pi, label='Algebraic Control Law')
ax_alpha_com.plot(dif_set[-1].t, dif_set[-1].u[0, :] * 180 / np.pi, '--', label='Differential Control Law')
# ax_alpha_com.legend()

fig.colorbar(plt.cm.ScalarMappable(norm=norm), cmap=cmap,
                cax=ax_cb, label=r'Initial FPA $\gamma_0$, [deg]',
                orientation='horizontal')

ax_path_alg.set_title('Algebraic Control Law: Path')
ax_path_alg.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_path_alg.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

# ax_hv_alg.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
# ax_hv_alg.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

ax_alpha_alg.set_title('Algebraic Control Law: AOA')
ax_alpha_alg.set_xlabel(r'Time $t$ [s]')
ax_alpha_alg.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')
# ax_alpha_alg.set_xlim(ax_alpha_dif.get_xlim())
# ax_alpha_alg.set_ylim(ax_alpha_dif.get_ylim())

ax_path_dif.set_title('Differential Control Law: Path')
ax_path_dif.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_path_dif.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')
ax_path_dif.set_xlim(ax_path_alg.get_xlim())
ax_path_dif.set_ylim(ax_path_alg.get_ylim())

# ax_hv_dif.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
# ax_hv_dif.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

ax_alpha_dif.set_title('Differential Control Law: AOA')
ax_alpha_dif.set_xlabel(r'Time $t$ [s]')
ax_alpha_dif.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

ax_path_com.set_title('Comparison: Path')
ax_path_com.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_path_com.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')
ax_path_com.set_xlim(ax_path_alg.get_xlim())
ax_path_com.set_ylim(ax_path_alg.get_ylim())

ax_alpha_com.set_title('Comparison: AOA')
ax_alpha_com.set_xlabel(r'Time $t$ [s]')
ax_alpha_com.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

fig.subplots_adjust(
        top=0.914,
        bottom=0.075,
        left=0.091,
        right=0.977,
        hspace=0.763,
        wspace=0.251
)

print(max(gam_0s_dif))

plt.show()
