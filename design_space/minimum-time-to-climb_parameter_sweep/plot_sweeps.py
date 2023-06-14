import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize
import scipy

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

isp_sweep = gp.load_sol_set('isp_series.data')
thrust_sweep = gp.load_sol_set('thrust_series.data')

fig = plt.figure(figsize=(6.5, 7.5))
title = fig.suptitle('Effect of Parameter Sweep\non Minimum Time-to-Climb Problem')

gs = gridspec.GridSpec(4, 2, height_ratios=(1, 1, 1, 0.05))

ax_isp_sweep = fig.add_subplot(gs[0, 0])
ax_thrust_sweep = fig.add_subplot(gs[0, 1])

ax_isp_traj = fig.add_subplot(gs[1, 0])
ax_thrust_traj = fig.add_subplot(gs[1, 1])

ax_isp_adj = fig.add_subplot(gs[2, 0])
ax_thrust_adj = fig.add_subplot(gs[2, 1])

cax_isp = fig.add_subplot(gs[3, 0])
cax_thrust = fig.add_subplot(gs[3, 1])

isp_list = []
tf_isp = []
int_lambda_w = []
isp_predictive_difference = []
for sol in isp_sweep:
    isp_list.append(sol.k[-1])
    tf_isp.append(sol.t[-1])
    int_lambda_w.append(scipy.integrate.simpson(sol.lam[3, :], sol.t))

ax_isp_sweep.plot(isp_list, tf_isp, linewidth=2)
ax_isp_sweep.set_xlabel(r'Specific Impulse, $I_{sp}$ [sec]')
ax_isp_sweep.set_ylabel(r'Final Time $t_f$ [sec]')

thrust_mult_list = []
tf_thrust = []
for sol in thrust_sweep:
    thrust_mult_list.append(sol.k[-1])
    tf_thrust.append(sol.t[-1])

ax_thrust_sweep.plot(thrust_mult_list, tf_thrust, linewidth=2)
ax_thrust_sweep.set_xlabel(r'Thrust Mult., $k_{\text{Thrust}}$')
ax_thrust_sweep.set_ylabel(r'Final Time $t_f$ [sec]')
# ax_trade_off.legend()

cmap = plt.colormaps['viridis_r']
cmap_r = plt.colormaps['viridis_r']

SKIP = 1

max_isp = isp_list[-1]
min_isp = isp_list[0]
for sol in isp_sweep[::]:
    ax_isp_traj.plot(sol.t, sol.x[0, :] / 1e3, color=cmap((sol.k[-1] - min_isp) / (max_isp - min_isp)), linewidth=1)
    ax_isp_adj.plot(sol.t, np.rad2deg(sol.u[0]), color=cmap((sol.k[-1] - min_isp) / (max_isp - min_isp)), linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_isp, max_isp), cmap=cmap),
             cax=cax_isp, label=r'Specific Impulse, $I_{sp}$ [sec]', orientation='horizontal', aspect=50)

# ax_isp_traj.set_title('Max Terminal Downrange')
ax_isp_traj.set_xlabel(r'Time, $t$ [s]')
ax_isp_traj.set_ylabel(r'Altitude, $h$ [$\times 10^3$ ft]')

ax_isp_adj.set_xlabel(r'Time, $t$ [s]')
ax_isp_adj.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

max_thrust_mult = thrust_mult_list[-1]
min_thrust_mult = thrust_mult_list[0]
for sol in thrust_sweep[::]:
    ax_thrust_traj.plot(sol.t, sol.x[0, :] / 1e3,
                        color=cmap((sol.k[-1] - min_thrust_mult) / (max_thrust_mult - min_thrust_mult)), linewidth=1)
    ax_thrust_adj.plot(sol.t, np.rad2deg(sol.u[0]),
                       color=cmap((sol.k[-1] - min_thrust_mult) / (max_thrust_mult - min_thrust_mult)), linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_thrust_mult, max_thrust_mult), cmap=cmap),
             cax=cax_thrust, label=r'Thrust Mult., $k_{\text{Thrust}}$', orientation='horizontal', aspect=50)

ax_thrust_traj.set_xlabel(r'Time, $t$ [s]')
ax_thrust_traj.set_ylabel(r'Altitude, $h$ [$\times 10^3$ ft]')

ax_thrust_adj.set_xlabel(r'Time, $t$ [s]')
ax_thrust_adj.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

fig.subplots_adjust(
        top=0.911,
        bottom=0.084,
        left=0.116,
        right=0.96,
        hspace=0.584,
        wspace=0.354
)

plt.show()
