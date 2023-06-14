import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec

rcParams['font.family'] = 'serif'
# rcParams['font.size'] = 10
rcParams['mathtext.fontset'] = 'stix'

crossrange_set = gp.load_sol_set('crossrange.data')
mitigated_set = gp.load_sol_set('crossrange_mit.data')

sol = crossrange_set[-1]
mit_sol = mitigated_set[-1]
# heat_sol = heat_rate_set[-1]

fig = plt.figure(figsize=(6.5, 7))
title = fig.suptitle('Space Shuttle Crossrange Problem')

gs = gridspec.GridSpec(3, 2)

ax_ground = fig.add_subplot(gs[0, 0])
# ax_heat = fig.add_subplot(gs[0, 1])
ax_alt = fig.add_subplot(gs[0, 1])

ax_hv = fig.add_subplot(gs[1, 0])
ax_gam = fig.add_subplot(gs[1, 1])

ax_alpha = fig.add_subplot(gs[2, 0])
ax_beta = fig.add_subplot(gs[2, 1])

fig_adj = plt.figure(figsize=(6.5, 7))

ax_lam_h = fig_adj.add_subplot(gs[0, 0])
ax_lam_phi = fig_adj.add_subplot(gs[0, 1])
ax_lam_theta = fig_adj.add_subplot(gs[1, 0])
ax_lam_v = fig_adj.add_subplot(gs[1, 1])
ax_lam_gam = fig_adj.add_subplot(gs[2, 0])
ax_lam_psi = fig_adj.add_subplot(gs[2, 1])


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


ax_ground.plot(sol.x[1, :] * 180 / np.pi, sol.x[2, :] * 180 / np.pi,
               label='No Path Constraint')
ax_ground.plot(mit_sol.x[1, :] * 180 / np.pi, mit_sol.x[2, :] * 180 / np.pi,
               label=r'$|\beta|$ < 70 deg')
# ax_ground.plot(heat_sol.x[1, :] * 180 / np.pi, heat_sol.x[2, :] * 180 / np.pi,
#                label=r'$q$ < 70 $\mathregular{BTU/ft^2/sec}$')
ax_ground.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_ground.set_ylabel(r'Crossrange, $\theta$ [deg]')

fig.legend(loc='lower center', ncols=3)

ax_hv.plot(sol.x[3, :] / 1e3, sol.x[0, :] / 1e5, label='No Path Constraint')
ax_hv.plot(mit_sol.x[3, :] / 1e3, mit_sol.x[0, :] / 1e5, label=r'$\left|\beta\right|$ < 70 deg')
# ax_hv.plot(heat_sol.x[3, :] / 1e3, heat_sol.x[0, :] / 1e5, label=r'$q < 70 \mathregular{BTU/ft^2/sec}$')
ax_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
ax_hv.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

# ax_alt.plot(compute_range(sol.x[1, :], sol.x[2, :]) * 20_902_900 / 5280, sol.x[0, :] / 1e5,
#             label='No Path Constraint')
# ax_alt.plot(compute_range(heat_sol.x[1, :], heat_sol.x[2, :]) * 20_902_900 / 5280, heat_sol.x[0, :] / 1e5,
#             label='q < 70 BTU/ft^2/sec')

ax_alt.plot(sol.t, sol.x[0, :] / 1e5, label='No Path Constraint')
ax_alt.plot(mit_sol.t, mit_sol.x[0, :] / 1e5, label='No Path Constraint')
# ax_alt.plot(heat_sol.t, heat_sol.x[0, :] / 1e5, label=r'$q < 70 \mathregular{BTU/ft^2/sec}$')
ax_alt.set_xlabel(r'Time $t$ [s]')
ax_alt.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

ax_gam.plot(sol.t, sol.x[4, :] * 180 / np.pi, label='No Path Constraint')
ax_gam.plot(mit_sol.t, mit_sol.x[4, :] * 180 / np.pi, label=r'$q < 70 \mathregular{BTU/ft^2/sec}$')
# ax_gam.plot(heat_sol.t, heat_sol.x[4, :] * 180 / np.pi, label=r'$q < 70 \mathregular{BTU/ft^2/sec}$')
ax_gam.set_xlabel(r'Time $t$ [s]')
ax_gam.set_ylabel(r'FPA, $\gamma$ [deg]')

# ax_heat.plot(sol.t, compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :]), label='No Path Constraint')
# ax_heat.plot(bank_sol.t, compute_heating(bank_sol.x[0, :], bank_sol.x[3, :], bank_sol.u[0, :]),
#              label=r'$|\beta|$ < 70 deg')
# ax_heat.plot(heat_sol.t, compute_heating(heat_sol.x[0, :], heat_sol.x[3, :], heat_sol.u[0, :]),
#              label=r'$q < 70 \mathregular{BTU/ft^2/sec}$')
# ax_heat.set_xlabel(r'Time $t$ [s]')
# ax_heat.set_ylabel(r'Heat-Rate,' + '\n' + r'$q$ [$\mathregular{BTU/ft^2/sec}$]')

ax_alpha.plot(sol.t, sol.u[0, :] * 180 / np.pi, '*-', label='No Path Constraint')
ax_alpha.plot(mit_sol.t, mit_sol.u[0, :] * 180 / np.pi, '+--', label=r'$|\beta|$ < 70 deg')
# ax_alpha.plot(heat_sol.t, heat_sol.u[0, :] * 180 / np.pi, label=r'$q < 70 \mathregular{BTU/ft^2/sec}$')
# ax_alpha.set_ylim(0, 30)
ax_alpha.set_xlabel(r'Time $t$ [s]')
ax_alpha.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

ax_beta.plot(sol.t, sol.u[1, :] * 180 / np.pi, label='No Path Constraint')
ax_beta.plot(mit_sol.t, mit_sol.u[1, :] * 180 / np.pi, label=r'$|\beta|$ < 70 deg')
# ax_beta.plot(heat_sol.t, heat_sol.u[1, :] * 180 / np.pi, label=r'$q < 70 \mathregular{BTU/ft^2/sec}$')
ax_beta.set_xlabel(r'Time $t$ [s]')
ax_beta.set_ylabel(r'Bank Angle, $\beta$ [deg]')

ax_lam_h.plot(sol.t, sol.lam[0, :])
ax_lam_phi.plot(sol.t, sol.lam[1, :])
ax_lam_theta.plot(sol.t, sol.lam[2, :])
ax_lam_v.plot(sol.t, sol.lam[3, :])
ax_lam_gam.plot(sol.t, sol.lam[4, :])
ax_lam_psi.plot(sol.t, sol.lam[5, :])

ax_lam_h.plot(mit_sol.t, mit_sol.lam[0, :])
ax_lam_phi.plot(mit_sol.t, mit_sol.lam[1, :])
ax_lam_theta.plot(mit_sol.t, mit_sol.lam[2, :])
ax_lam_v.plot(mit_sol.t, mit_sol.lam[3, :])
ax_lam_gam.plot(mit_sol.t, mit_sol.lam[4, :])
ax_lam_psi.plot(mit_sol.t, mit_sol.lam[5, :])

fig.subplots_adjust(
        top=0.925,
        bottom=0.125,
        left=0.125,
        right=0.975,
        hspace=0.4,
        wspace=0.4
)

plt.show()
