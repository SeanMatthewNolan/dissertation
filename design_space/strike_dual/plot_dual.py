import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec,rc
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

max_phi = gp.load_sol_set('max_phi.bin')
max_v = gp.load_sol_set('max_v.bin')

fig = plt.figure(figsize=(6.5, 7.5))
title = fig.suptitle('Equivalence of Maximum Downrange\nand Maximum Terminal Velocity Problems')

gs = gridspec.GridSpec(5, 1, height_ratios=(1, 1, 0.1, 1, 0.1))

ax_trade_off = fig.add_subplot(gs[0, 0])
ax_max_phi = fig.add_subplot(gs[1, 0])
ax_cb_phi = fig.add_subplot(gs[2, 0])
ax_max_v = fig.add_subplot(gs[3, 0])
ax_cb_v = fig.add_subplot(gs[4, 0])

max_phi_phi = []
max_phi_v = []
for sol in max_phi:
    max_phi_phi.append(sol.x[1, -1])
    max_phi_v.append(sol.x[2, -1])

max_v_phi = []
max_v_v = []
for sol in max_v:
    max_v_phi.append(sol.x[1, -1])
    max_v_v.append(sol.x[2, -1])

ax_trade_off.plot(np.rad2deg(max_phi_phi), np.array(max_phi_v) / 1e3, label=r'Max $\phi_f$', linewidth=2)
ax_trade_off.plot(np.rad2deg(max_v_phi), np.array(max_v_v) / 1e3, ':', label=r'Max $v_f$', linewidth=3)
ax_trade_off.set_xlabel(r'Terminal Downrange Angle, $\phi_f$ [deg]')
ax_trade_off.set_ylabel('Terminal Velocity,\n' + r'$v_f$ [$\times10^3$ ft/s]')
ax_trade_off.legend()


cmap = plt.colormaps['viridis']
cmap_v = plt.colormaps['viridis_r']

SKIP = 1

max_vf = max_phi[-1].x[2, -1]
min_vf = max_phi[0].x[2, -1]
for sol in max_phi[::-SKIP]:
    ax_max_phi.plot(np.rad2deg(sol.x[1, :]), sol.x[0, :] / 1e5,
                    color=cmap((sol.x[2, -1] - min_vf) / (max_vf - min_vf)), linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_vf, max_vf), cmap=cmap),
             cax=ax_cb_phi, label='Constrained Terminal Velocity [ft/s]', orientation='horizontal', aspect=50)

ax_max_phi.set_title('Max Terminal Downrange')
ax_max_phi.set_xlabel('Downrange Angle, $\phi$ [deg]')
ax_max_phi.set_ylabel('Altitude,\n' + r'$h$ [$\times 10^5$ ft]')

max_phif = np.rad2deg(max_v[-1].x[1, -1])
min_phif = np.rad2deg(max_v[0].x[1, -1])
for sol in max_v[::SKIP]:
    ax_max_v.plot(np.rad2deg(sol.x[1, :]), sol.x[0, :] / 1e5,
                    color=cmap_v((np.rad2deg(sol.x[1, -1]) - min_phif) / (max_phif - min_phif)), linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_phif, max_phif), cmap=cmap_v),
             cax=ax_cb_v, label='Constrained Terminal Downrange [deg]', orientation='horizontal', aspect=50)

ax_max_v.set_title('Max Terminal Velocity')
ax_max_v.set_xlabel('Downrange Angle, $\phi$ [deg]')
ax_max_v.set_ylabel('Altitude,\n' + r'$h$ [$\times 10^5$ ft]')

fig.subplots_adjust(
        top=0.917,
        bottom=0.073,
        left=0.129,
        right=0.953,
        hspace=1.0,
        wspace=0.2
)

plt.show()
