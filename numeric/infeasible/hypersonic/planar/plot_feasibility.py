import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# max_phi = gp.load_sol_set('max_phi.bin')
max_v = gp.load_sol_set('max_v.bin')
max_v_r = gp.load_sol_set('max_v_reverse.bin')

max_v.solutions = max_v_r.solutions[::-1] + max_v.solutions

fig = plt.figure(figsize=(6.5, 7.5))
title = fig.suptitle('Adjoint Variables and Feasibility Boundary:\nMaximum Terminal Velocity Problems')

gs = gridspec.GridSpec(4, 1, height_ratios=(1, 1, 1, 0.1))

ax_adjoint = fig.add_subplot(gs[0, 0])
ax_trade_off = fig.add_subplot(gs[1, 0])
ax_max_v = fig.add_subplot(gs[2, 0])
ax_cb_v = fig.add_subplot(gs[3, 0])

max_v_phi = []
max_v_v = []
nus = []
adj_norm = []
for sol in max_v:
    max_v_phi.append(sol.x[1, -1])
    max_v_v.append(sol.x[2, -1])
    nus.append(np.concatenate((sol.nu0, sol.nuf)))
    adj_norm.append(np.linalg.norm(np.concatenate((sol.nu0, sol.nuf))))

nus = np.array(nus)

# ax_trade_off.plot(np.rad2deg(max_v_phi), nus, label=r'Max $v_f$')
ax_trade_off.plot(np.rad2deg(max_v_phi), np.array(max_v_v) / 1e3, label=r'Max $v_f$')
ax_trade_off.set_xlabel(r'Terminal Downrange Angle, $\phi_f$ [deg]')
ax_trade_off.set_ylabel('Terminal Velocity,\n' + r'$v_f$ [$\times10^3$ ft/s]')
# ax_trade_off.legend()


cmap = plt.colormaps['viridis']
cmap_v = plt.colormaps['viridis_r']
max_phif = np.rad2deg(max_v[-1].x[1, -1])
min_phif = np.rad2deg(max_v[0].x[1, -1])

ax_adjoint.plot(np.rad2deg(max_v_phi), np.array(adj_norm))

ax_adjoint.set_title('Norm of Static Adjoint Variables')
# ax_adjoint.set_yscale('log')
ax_adjoint.set_xlabel('Terminal Downrange Angle, $\phi$ [deg]')
ax_adjoint.set_ylabel(r'Static Adjoints Norm, $|| \boldsymbol{\nu} ||$')

for sol in max_v[::]:

    ax_max_v.plot(np.rad2deg(sol.x[1, :]), sol.x[0, :] / 1e5,
                  color=cmap_v((np.rad2deg(sol.x[1, -1]) - min_phif) / (max_phif - min_phif)), linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_phif, max_phif), cmap=cmap_v),
             cax=ax_cb_v, label='Constrained Terminal Downrange [deg]', orientation='horizontal', aspect=50)

ax_max_v.set_title('Max Terminal Velocity')
ax_max_v.set_xlabel('Downrange Angle, $\phi$ [deg]')
ax_max_v.set_ylabel('Altitude,\n' + r'$h$ [$\times 10^5$ ft]')

fig.subplots_adjust(
        top=0.891,
        bottom=0.075,
        left=0.132,
        right=0.977,
        hspace=0.805,
        wspace=0.2
)

plt.show()
