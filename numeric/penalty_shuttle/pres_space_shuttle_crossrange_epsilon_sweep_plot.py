import giuseppe as gp
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc
from matplotlib.colors import Normalize, LogNorm

import scipy


rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

eps_sweep = gp.load_sol_set('eps_effect.data')
# eps_sweep = gp.load_sol_set('eps_effect_sec.data')
# eps_sweep = gp.load_sol_set('eps_effect_slight.data')
# eps_sweep = gp.load_sol_set('eps_effect_inactive.data')

fig = plt.figure(figsize=(3.5, 2.75))
title = fig.suptitle('Effects of Smoothing Constant\nSpace Shuttle Crossrange Problem')

gs = gridspec.GridSpec(2, 2, height_ratios=(1, 0.2))

ax_cost = fig.add_subplot(gs[0, 0])
ax_rel_cost_sens = fig.add_subplot(gs[0, 1])


def compute_density(altitude):
    return 0.002378 * np.exp(-altitude / 23800)


def compute_heating(altitude, velocity, alpha):
    alpha_hat = alpha * 180 / np.pi
    rho = compute_density(altitude)
    q_r = 17700 * np.sqrt(rho) * (0.0001 * velocity) ** 3.07
    q_a = 1.0672181 + -0.19213774e-1 * alpha_hat + 0.21286289e-3 * alpha_hat**2 + -0.10117249e-5 * alpha_hat**3
    return q_r * q_a


cmap = plt.colormaps['viridis_r']

buffers = []
pred_diff = []
int_mu = []
max_dp_ds = []
min_mu = []
mean_mu = []
max_mu = []
epsilons = []
cross_f = []
costs = []
penalties = []
for k, sol in enumerate(eps_sweep):
    q = compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :])
    cross_f.append(-(sol.x[2, -1]))
    costs.append(sol.cost)
    epsilon = sol.k[25]
    q_u = sol.k[26]
    epsilons.append(epsilon)
    int_mu.append(scipy.integrate.simpson(sol.aux['mu Penalty: q'], sol.t))
    min_mu.append(np.min(np.abs(sol.aux['mu Penalty: q'])))
    mean_mu.append(np.mean(np.abs(sol.aux['mu Penalty: q'])))
    max_mu.append(np.max(np.abs(sol.aux['mu Penalty: q'])))
    penalties.append(scipy.integrate.simpson(sol.aux['Penalty: q'], sol.t))
    buffers.append(min(q_u - q))

int_mu = np.asarray(int_mu)
min_mu = np.asarray(min_mu)
mean_mu = np.asarray(mean_mu)
max_mu = np.asarray(max_mu)
costs = np.asarray(costs)
cross_f = np.asarray(cross_f)
epsilons = np.asarray(epsilons)
penalties = np.asarray(penalties)
buffers = np.asarray(buffers)

for sol in eps_sweep[::]:
    color = cmap((np.log10(sol.k[25]) - np.log10(min(epsilons))) / (np.log10(max(epsilons)) - np.log10(min(epsilons))))

    epsilon = sol.k[25]
    q_max = sol.k[26]
    q = compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :])
    dq_dt = np.diff(q) / np.diff(sol.t)


ax_cost.plot(epsilons, costs, label=r'Reg. Cost $J$', linewidth=2)
ax_cost.plot(epsilons, cross_f, label=r'Orig. Cost $-\theta_f$', linewidth=2)
ax_cost.set_xscale('log')
ax_cost.set_xlabel(r'$\epsilon_q$ [rad]')
ax_cost.set_ylabel(r'Cost, $\bar{J}$ [rad]')


ax_rel_cost_sens.plot(epsilons, max_mu / epsilons, label=r'Max.', linewidth=2)
ax_rel_cost_sens.plot(epsilons, mean_mu / epsilons, label=r'Mean', linewidth=2)

ax_rel_cost_sens.set_xscale('log')
ax_rel_cost_sens.set_yscale('log')
ax_rel_cost_sens.set_xlabel(r'$\epsilon_q$ [rad]')
ax_rel_cost_sens.set_ylabel(
        'Penalty Sensitivity \n'
        + r'$\frac{\partial P}{\partial C}$,'
        + r'[$\mathrm{{ft}^2 / {BTU}}$]')

fig.subplots_adjust(
        top=0.794,
        bottom=0.235,
        left=0.198,
        right=0.939,
        hspace=0.691,
        wspace=0.982
)

fig.legend(*ax_cost.get_legend_handles_labels(), bbox_to_anchor=gs[1, 0].get_position(fig), loc='upper center')
fig.legend(*ax_rel_cost_sens.get_legend_handles_labels(), bbox_to_anchor=gs[1, 1].get_position(fig), loc='upper center')

plt.show()
