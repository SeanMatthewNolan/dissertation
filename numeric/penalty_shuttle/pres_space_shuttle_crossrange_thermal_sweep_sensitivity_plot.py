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

fig = plt.figure(figsize=(3.5, 2.75))
title = fig.suptitle('Cost Sensitivity to Maximum Heat Rate\nSpace Shuttle Crossrange Problem')

gs = gridspec.GridSpec(2, 2, height_ratios=(1, 0.2))

ax_cost = fig.add_subplot(gs[0, 0])
ax_cost_sens = fig.add_subplot(gs[0, 1])

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

ax_cost.plot(q_us, cross_f, label=r'Original Cost $-\theta_f$', linewidth=2, color='C0')
ax_cost.plot(q_us, cost, '--', label=r'Augmented Cost $J$', linewidth=2, color='C1')
ax_cost.set_xlabel(r'$q_U$, [$\mathrm{BTU/ft^2/sec}$]')
ax_cost.set_ylabel(r'Cost $J$ [deg]')

ax_cost_sens.plot(
        q_us, np.rad2deg(pred_diff), label=r'$\int_0^{t_f} \epsilon \frac{\partial P}{\partial C} \, dt$',
        linewidth=2
)
ax_cost_sens.plot(
        (np.asarray(q_us[:-1]) + np.asarray(q_us[1:])) / 2, np.diff(cross_f) / np.diff(q_us), '--',
        linewidth=2, label=r'$\frac{\Delta J}{\Delta q_U}$'
)
ax_cost_sens.set_xlabel(r'$q_U$, [$\mathrm{BTU/ft^2/sec}$]')
ax_cost_sens.set_ylabel(r'Cost Sens. $-\frac{\partial\theta_f}{\partial q_U}$')

fig.subplots_adjust(
        top=0.794,
        bottom=0.235,
        left=0.198,
        right=0.939,
        hspace=0.691,
        wspace=0.982
)

fig.legend(*ax_cost.get_legend_handles_labels(), bbox_to_anchor=gs[1, 0].get_position(fig), loc='upper center')
fig.legend(*ax_cost_sens.get_legend_handles_labels(), bbox_to_anchor=gs[1, 1].get_position(fig), loc='upper center')

plt.show()
