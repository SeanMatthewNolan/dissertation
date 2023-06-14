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

fig = plt.figure(figsize=(6.5, 7.5))
title = fig.suptitle('Effects of Smoothing Constant\nSpace Shuttle Crossrange Problem')

gs = gridspec.GridSpec(5, 2, height_ratios=(1, 1, 1, 1, 0.1))

ax_heat = fig.add_subplot(gs[0, 0])
ax_cost = fig.add_subplot(gs[0, 1])
ax_mu = fig.add_subplot(gs[1, 0])
ax_cost_sens = fig.add_subplot(gs[1, 1])
ax_penalty = fig.add_subplot(gs[2, 0])
ax_rel_cost_sens = fig.add_subplot(gs[2, 1])
ax_buffer = fig.add_subplot(gs[3, 0])
ax_min_buffer = fig.add_subplot(gs[3, 1])

cax = fig.add_subplot(gs[-1, :])


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

    ax_heat.plot(sol.t, q, color=color)
    ax_mu.plot(sol.t, sol.aux['mu Penalty: q'] , color=color)
    ax_penalty.plot(sol.t, sol.aux['mu Penalty: q'] / epsilon, color=color)
    # ax_penalty.plot(sol.t[1:], sol.aux['mu Penalty: q'][1:] * dq_dt / epsilon, color=color)
    ax_buffer.plot(sol.t, (q_max - q), color=color)


ax_heat.set_xlabel(r'Time $t$ [s]')
ax_heat.set_ylabel('Heat-Rate, $q$,\n' + r'[$\mathrm{BTU/ft^2/sec}$]')

ax_mu.set_xlabel(r'Time $t$ [s]')
ax_mu.set_ylabel(
        'Approx. Adjoint\n'
        + r'$\tilde{\mu} = \epsilon \frac{\partial P}{\partial C}$,'
        + '\n' + r'[$\mathrm{{rad} \; {ft}^2 / {BTU}}$]'
)

ax_penalty.set_xlabel(r'Time $t$ [s]')
ax_penalty.set_ylabel(
        'Penalty Sensitivity\n'
        + r'$\frac{\partial P}{\partial C}$,'
        + r'[$\mathrm{{ft}^2 / {BTU}}$]')

ax_buffer.set_yscale('log')
ax_buffer.set_xlabel(r'Time $t$ [s]')
ax_buffer.set_ylabel('Dist. from Boundary,\n' + r'$q_U - q$' + '\n' + '[$\mathrm{BTU/ft^2/sec}$]')

ax_cost.plot(epsilons, costs, label=r'Reg. Cost $\bar{J}$', linewidth=2)
ax_cost.plot(epsilons, cross_f, label=r'Orig. Cost $-\theta_f$', linewidth=2)
ax_cost.set_xscale('log')
ax_cost.set_xlabel(r'Smoothing Constant, $\epsilon_q$ [rad]')
ax_cost.set_ylabel(r'Cost, $\bar{J}$ [rad]')
ax_cost.legend()

ax_cost_sens.plot(epsilons, np.asarray(int_mu), linewidth=2)
# ax_cost_sens.plot(epsilons, np.asarray(int_mu), linewidth=2)

# ax_cost_sens.legend()
ax_cost_sens.set_xscale('log')
ax_cost_sens.set_xlabel(r'Smoothing Constant, $\epsilon_q$ [rad]')
ax_cost_sens.set_ylabel('Cost Sensitivity,\n'
                        + r'$\frac{\partial J}{\partial q}$ [$\mathrm{rad \; s \; ft^2 /BTU}$]')

# ax_rel_cost_sens.plot(epsilons, max_mu / epsilons, label=r'$\max{\frac{\partial P}{\partial q}}$', linewidth=2)
# ax_rel_cost_sens.plot(epsilons, mean_mu / epsilons,
#                       label=r'$\left(\left| \int_0^{t_f} \frac{\partial P}{\partial q} \; dt \right|\right) / t_f$',
#                       linewidth=2)
ax_rel_cost_sens.plot(epsilons, max_mu / epsilons, label=r'Max.', linewidth=2)
ax_rel_cost_sens.plot(epsilons, mean_mu / epsilons, label=r'Mean', linewidth=2)

ax_rel_cost_sens.legend()
ax_rel_cost_sens.set_xscale('log')
ax_rel_cost_sens.set_yscale('log')
ax_rel_cost_sens.set_xlabel(r'Smoothing Constant, $\epsilon_q$ [rad]')
ax_rel_cost_sens.set_ylabel(
        'Penalty Sensitivity \n'
        + r'$\frac{\partial P}{\partial C}$,'
        + r'[$\mathrm{{ft}^2 / {BTU}}$]')
print(linregress(np.log(epsilons), np.log(max_mu / epsilons)))

ax_min_buffer.plot(epsilons, buffers)
ax_min_buffer.set_xlabel(r'Smoothing Constant, $\epsilon_q$ [rad]')
ax_min_buffer.set_ylabel(r'Min. Dist.' + '\n' + r'$\min{(q_U - q)}$' + '\n' + '[$\mathrm{BTU/ft^2/sec}$]')
ax_min_buffer.set_xscale('log')
ax_min_buffer.set_yscale('log')

print(linregress(np.log(epsilons), np.log(buffers)))

fig.colorbar(plt.cm.ScalarMappable(norm=LogNorm(min(epsilons), max(epsilons)), cmap=cmap),
             cax=cax, label=r'Smoothing Constant, $\epsilon_q$ [rad]',
             orientation='horizontal', aspect=50)

fig.tight_layout()

# fig.subplots_adjust(
#         top=0.918,
#         bottom=0.076,
#         left=0.166,
#         right=0.977,
#         hspace=0.53,
#         wspace=0.302
# )

plt.show()
