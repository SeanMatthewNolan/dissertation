import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc
from matplotlib.colors import Normalize

import scipy

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

fpa_sweep = gp.load_sol_set('heat_rate.data')

fig = plt.figure(figsize=(6.5, 7.5))
title = fig.suptitle('Cost Sensitivity to Maximum Heat Rate\nSpace Shuttle Crossrange Problem')

gs = gridspec.GridSpec(4, 2, height_ratios=(1, 1, 1, 0.05))

ax_cost = fig.add_subplot(gs[1, 1])
# ax_cost_diff = fig.add_subplot(gs[1:3, :])
ax_cost_sens = fig.add_subplot(gs[0, :])
ax_cost_sens_zoom = fig.add_subplot(gs[1, 0])

# ax_hv = fig.add_subplot(gs[1, 0])
ax_heat = fig.add_subplot(gs[2, 0])
ax_ground = fig.add_subplot(gs[2, 1])
# ax_alt = fig.add_subplot(gs[1, 0])

# ax_alpha = fig.add_subplot(gs[2, 0])
# ax_beta = fig.add_subplot(gs[2, 1])

cax = fig.add_subplot(gs[3, :])


def haversine(_ang):
    return (1 - np.cos(_ang)) / 2


def compute_range(downrange, crossrange):
    return np.arctan(haversine(crossrange) + np.cos(crossrange) * haversine(downrange))


def compute_density(altitude):
    return 0.002378 * np.exp(-altitude / 23800)


def compute_heating(altitude, velocity, alpha):
    alpha_hat = alpha * 180 / np.pi
    rho = compute_density(altitude)
    q_r = 17700 * np.sqrt(rho) * (0.0001 * velocity) ** 3.07
    q_a = 1.0672181 + -0.19213774e-1 * alpha_hat + 0.21286289e-3 * alpha_hat**2 + -0.10117249e-5 * alpha_hat**3
    return q_r * q_a


cmap = plt.colormaps['viridis_r']


def p(_s):
    return 1 / _s


def dp_ds(_s):
    return -1 / _s**2


pred_diff = []
q_us = []
cross_f = []
cost = []
for sol in fpa_sweep:
    q = compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :])
    cross_f.append(-np.rad2deg(sol.x[2, -1]))
    cost.append(np.rad2deg(sol.cost))
    epsilon = sol.k[25]
    q_u = sol.k[26]
    q_us.append(q_u)
    s = q_u - q
    pred_diff.append(
            scipy.integrate.simpson(epsilon * dp_ds(s), sol.t)
    )

for sol in fpa_sweep[::-5]:
    color = cmap((max(compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :])) - min(q_us)) / (max(q_us) - min(q_us)))

    ax_ground.plot(sol.x[1, :] * 180 / np.pi, sol.x[2, :] * 180 / np.pi, color=color)

    # ax_hv.plot(sol.x[3, :] / 1e3, sol.x[0, :] / 1e5, color=color)
    # ax_heat.plot(sol.t, sol.x[4, :] * 180 / np.pi, color=color)

    ax_heat.plot(sol.t, compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :]), color=color)

    # ax_alpha.plot(sol.t, sol.u[0, :] * 180 / np.pi, color=color)
    # ax_beta.plot(sol.t, sol.u[1, :] * 180 / np.pi, color=color)


ax_ground.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_ground.set_ylabel(r'Crossrange, $\theta$ [deg]')

# ax_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
# ax_hv.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

ax_heat.set_xlabel(r'Time $t$ [s]')
ax_heat.set_ylabel('Heat-Rate, $q$,\n' + r'[$\mathrm{BTU/ft^2/sec}$]')

ax_cost.plot(q_us, cross_f, ':', label=r'Original Cost $\bar{J}$', linewidth=2, color='C2')
ax_cost.plot(q_us, cost, label=r'Augmented Cost $J$', linewidth=2, color='C0')
ax_cost.legend()
ax_cost.set_xlabel(r'Maximum Heat-Rate, $q_U$, [$\mathrm{BTU/ft^2/sec}$]')
ax_cost.set_ylabel(r'Cost, $-\theta_f$ [deg]')

ax_cost_sens.plot(
        q_us, np.rad2deg(pred_diff), label=r'$\int_0^{t_f} \epsilon \frac{\partial P}{\partial S} \, dt$',
        linewidth=2
)
ax_cost_sens.plot(
        (np.asarray(q_us[:-1]) + np.asarray(q_us[1:])) / 2, np.diff(cost) / np.diff(q_us), '--',
        linewidth=2, label=r'$\frac{\Delta \bar{J}}{\Delta q_U}$'
)
ax_cost_sens.plot(
        (np.asarray(q_us[:-1]) + np.asarray(q_us[1:]))/2, np.diff(cross_f) / np.diff(q_us), ':',
        linewidth=2, label=r'$\frac{\Delta J}{\Delta q_U}$'
)
# ax_cost_sens.plot(
#         q_us, (cross_f - np.asarray(cost)), '-.',
#         linewidth=2, label='Penalty'
# )
ax_cost_sens.legend()
ax_cost_sens.set_xlabel(r'Maximum Heat-Rate, $q_U$, [$\mathrm{BTU/ft^2/sec}$]')
ax_cost_sens.set_ylabel('Cost Sensitivity,\n'
                        + r'$-\frac{\partial\theta_f}{\partial q_U}$ [$\mathrm{deg s ft^2 /BTU}$]')

zoom_slice = slice(-30, -1)
ax_cost_sens_zoom.plot(
        q_us[zoom_slice], np.rad2deg(pred_diff)[zoom_slice],
        label=r'$\int_0^{t_f} \epsilon \frac{\partial P}{\partial S} \, dt$',
        linewidth=2
)
ax_cost_sens_zoom.plot(
        ((np.asarray(q_us[:-1]) + np.asarray(q_us[1:])) / 2)[zoom_slice], (np.diff(cost) / np.diff(q_us))[zoom_slice],
        '--',
        linewidth=2, label=r'$\frac{\Delta \bar{J}}{\Delta q_U}$'
)
ax_cost_sens_zoom.plot(
        ((np.asarray(q_us[:-1]) + np.asarray(q_us[1:]))/2)[zoom_slice], (np.diff(cross_f) / np.diff(q_us))[zoom_slice],
        ':',
        linewidth=2, label=r'$\frac{\Delta J}{\Delta q_U}$'
)

ax_cost_sens_zoom.legend()
ax_cost_sens_zoom.set_xlabel(r'Maximum Heat-Rate, $q_U$, [$\mathrm{BTU/ft^2/sec}$]')
ax_cost_sens_zoom.set_ylabel('Cost Sensitivity,\n'
                             + r'$-\frac{\partial\theta_f}{\partial q_U}$ [$\mathrm{deg s ft^2 /BTU}$]')


# ax_alpha.set_xlabel(r'Time $t$ [s]')
# ax_alpha.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

# ax_beta.set_xlabel(r'Time $t$ [s]')
# ax_beta.set_ylabel(r'Bank Angle, $\beta$ [deg]')

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min(q_us), max(q_us)), cmap=cmap),
             cax=cax, label=r'Maximum Heat-Rate, $q_U$, [$\mathrm{BTU/ft^2/sec}$]',
             orientation='horizontal', aspect=50)

fig.subplots_adjust(
        top=0.918,
        bottom=0.076,
        left=0.166,
        right=0.977,
        hspace=0.53,
        wspace=0.302
)

plt.show()
