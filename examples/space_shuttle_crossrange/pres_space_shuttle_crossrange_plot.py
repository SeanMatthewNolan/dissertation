import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

crossrange_set = gp.load_sol_set('crossrange.data')
bank_set = gp.load_sol_set('crossrange_bank_constraint.data')
heat_rate_set = gp.load_sol_set('heat_rate.data')

sol = crossrange_set[-1]
bank_sol = bank_set[-1]
heat_sol = heat_rate_set[-1]

fig = plt.figure(figsize=(5.75, 2))
# title = fig.suptitle('Space Shuttle Crossrange Problem')

gs = gridspec.GridSpec(1, 4)

ax_ground = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[0, 1])

ax_alpha = fig.add_subplot(gs[0, 2])
ax_beta = fig.add_subplot(gs[0, 3])


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
ax_ground.plot(bank_sol.x[1, :] * 180 / np.pi, bank_sol.x[2, :] * 180 / np.pi,
               label=r'$\lvert\beta\lvert < 70$ deg')
ax_ground.plot(heat_sol.x[1, :] * 180 / np.pi, heat_sol.x[2, :] * 180 / np.pi,
               label=r'$q < 70 \text{BTU}/\text{ft}^2/\text{sec}$')
ax_ground.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_ground.set_ylabel(r'Crossrange, $\theta$ [deg]')

fig.legend(loc='lower center', ncols=3, frameon=False)

ax_heat.plot(sol.t, compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :]), label='No Path Constraint')
ax_heat.plot(bank_sol.t, compute_heating(bank_sol.x[0, :], bank_sol.x[3, :], bank_sol.u[0, :]),
             label=r'$|\beta|$ < 70 deg')
ax_heat.plot(heat_sol.t, compute_heating(heat_sol.x[0, :], heat_sol.x[3, :], heat_sol.u[0, :]),
             label=r'$q < 70 \text{BTU/ft^2/sec}$')
ax_heat.set_xlabel(r'Time $t$ [s]')
ax_heat.set_ylabel(r'Heat, $q$ [$\text{BTU}/\text{ft}^2/\text{s}$]')

ax_alpha.plot(sol.t, sol.u[0, :] * 180 / np.pi, label='No Path Constraint')
ax_alpha.plot(bank_sol.t, bank_sol.u[0, :] * 180 / np.pi, label=r'$|\beta|$ < 70 deg')
ax_alpha.plot(heat_sol.t, heat_sol.u[0, :] * 180 / np.pi, label=r'$q < 70 \text{BTU}/\text{ft}^2/\text{sec}$')
ax_alpha.set_xlabel(r'Time $t$ [s]')
ax_alpha.set_ylabel(r'AOA, $\alpha$ [deg]')

ax_beta.plot(sol.t, sol.u[1, :] * 180 / np.pi, label='No Path Constraint')
ax_beta.plot(bank_sol.t, np.arctan(bank_sol.u[1, :] / bank_sol.k[12]) * 2 / np.pi * bank_sol.k[14] * 180 / np.pi,
             label=r'$|\beta|$ < 70 deg')
ax_beta.plot(heat_sol.t, heat_sol.u[1, :] * 180 / np.pi, label=r'$q < 70 \text{BTU/ft^2/sec}$')
ax_beta.set_xlabel(r'Time $t$ [s]')
ax_beta.set_ylabel(r'Bank, $\beta$ [deg]')

fig.subplots_adjust(
        top=0.925,
        bottom=0.42,
        left=0.101,
        right=0.966,
        hspace=0.4,
        wspace=1.0
)

plt.show()
