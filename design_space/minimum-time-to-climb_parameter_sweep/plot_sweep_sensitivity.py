import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc
from matplotlib.colors import Normalize
import scipy

from giuseppe.utils.examples import Atmosphere1976

from lookup_tables import thrust_table_bspline, eta_table_bspline_expanded, CLalpha_table_bspline_expanded,\
    CD0_table_bspline_expanded, temp_table_bspline, dens_table_bspline


rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

isp_sweep = gp.load_sol_set('isp_series.data')
thrust_sweep = gp.load_sol_set('thrust_series.data')

fig = plt.figure(figsize=(6.5, 5))
title = fig.suptitle('Minimum Time-to-Climb\nDemonstration of Costate Relationship to Cost Functional')

gs = gridspec.GridSpec(3, 2, height_ratios=(1, 0.05, 1))

ax_isp_traj = fig.add_subplot(gs[0, 0])
ax_isp_control = fig.add_subplot(gs[0, 1])
cax_isp = fig.add_subplot(gs[1, :])
ax_isp_cost = fig.add_subplot(gs[2, 0])
ax_isp_diff = fig.add_subplot(gs[2, 1])

isp_list = []
tf_isp = []
int_lambda_w = []
isp_predictive_difference = []
for sol in isp_sweep:
    isp_list.append(sol.k[-1])
    tf_isp.append(sol.t[-1])
    int_lambda_w.append(scipy.integrate.simpson(sol.lam[3, :], sol.t))

    h = sol.x[0, :]
    v = sol.x[1, :]

    atm = Atmosphere1976(use_metric=False)

    # T = np.asarray([atm.temperature(alt) for alt in h])
    # rho = np.asarray([atm.density(alt) for alt in h])
    T = np.asarray(temp_table_bspline(h)).flatten()
    rho = np.asarray(dens_table_bspline(h)).flatten()

    a = np.sqrt(atm.specific_heat_ratio * atm.gas_constant * T)

    M = v / a
    Qdyn = 0.5 * rho * v ** 2

    thrust = np.asarray(thrust_table_bspline(np.vstack((M.T, h.T)))).flatten()

    isp_predictive_difference.append(
            scipy.integrate.simpson(sol.lam[3, :] * thrust / sol.k[-1]**2, sol.t)
    )

ax_isp_cost.plot(isp_list, tf_isp, linewidth=2)
ax_isp_cost.set_xlabel(r'Specific Impulse, $I_{sp}$ [sec]')
ax_isp_cost.set_ylabel(r'Final Time $t_f$ [sec]')

ax_isp_diff.plot((np.asarray(isp_list[:-1]) + np.asarray(isp_list[1:]))/2, np.diff(tf_isp) / np.diff(isp_list),
                 linewidth=2, label='From Continuation')
ax_isp_diff.plot(isp_list, isp_predictive_difference, ':', linewidth=2,
                 label=r'$\int_0^{t_f} \frac{T}{I_{sp}^2} \lambda_w \, dt$')
ax_isp_diff.set_xlabel(r'Specific Impulse, $I_{sp}$ [sec]')
ax_isp_diff.set_ylabel(r'Cost Sensitivity $\frac{d t_f}{d I_{sp}}$')
ax_isp_diff.legend()

cmap = plt.colormaps['viridis_r']

SKIP = 1

max_isp = isp_list[-1]
min_isp = isp_list[0]
for sol in isp_sweep[::]:
    ax_isp_traj.plot(sol.t, sol.x[0, :] / 1e3, color=cmap((sol.k[-1] - min_isp) / (max_isp - min_isp)), linewidth=1)
    ax_isp_control.plot(sol.t, np.rad2deg(sol.u[0, :]), color=cmap((sol.k[-1] - min_isp) / (max_isp - min_isp)),
                        linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_isp, max_isp), cmap=cmap),
             cax=cax_isp, label=r'Specific Impulse, $I_{sp}$ [sec]', orientation='horizontal', aspect=50)

# ax_isp_traj.set_title('Max Terminal Downrange')
ax_isp_traj.set_xlabel(r'Time, $t$ [s]')
ax_isp_traj.set_ylabel(r'Altitude, $h$ [$\times 10^3$ ft]')

ax_isp_control.set_xlabel(r'Time, $t$ [s]')
ax_isp_control.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')


fig.subplots_adjust(
        top=0.876,
        bottom=0.114,
        left=0.101,
        right=0.956,
        hspace=0.638,
        wspace=0.332
)

plt.show()
