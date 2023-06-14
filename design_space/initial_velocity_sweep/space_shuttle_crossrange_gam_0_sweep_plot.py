import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

fpa_sweep = gp.load_sol_set('initial_fpa_sweep.data')

fig = plt.figure(figsize=(6.5, 6.5))
title = fig.suptitle('Variation in Initial Flight-Path-Angle\nSpace Shuttle Crossrange Problem')

gs = gridspec.GridSpec(4, 2, height_ratios=(1, 1, 1, 0.05))

ax_cost = fig.add_subplot(gs[0, 0])
ax_ground = fig.add_subplot(gs[0, 1])

ax_hv = fig.add_subplot(gs[1, 0])
ax_gam = fig.add_subplot(gs[1, 1])
# ax_alt = fig.add_subplot(gs[1, 0])

ax_alpha = fig.add_subplot(gs[2, 0])
ax_beta = fig.add_subplot(gs[2, 1])

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

gamma_0 = []
cross_f = []
for sol in fpa_sweep:
    gamma_0.append(np.rad2deg(sol.x[4, 0]))
    cross_f.append(np.rad2deg(sol.x[2, -1]))

for sol in fpa_sweep[::-5]:
    color = cmap((np.rad2deg(sol.x[4, 0]) - min(gamma_0)) / (max(gamma_0) - min(gamma_0)))

    ax_ground.plot(sol.x[1, :] * 180 / np.pi, sol.x[2, :] * 180 / np.pi,
                   label='No Path Constraint', color=color)

    ax_hv.plot(sol.x[3, :] / 1e3, sol.x[0, :] / 1e5, label='No Path Constraint', color=color)

    ax_gam.plot(sol.t, sol.x[4, :] * 180 / np.pi, label='No Path Constraint', color=color)

    # ax_cost.plot(sol.t, compute_heating(sol.x[0, :], sol.x[3, :], sol.u[0, :]), label='No Path Constraint')

    ax_alpha.plot(sol.t, sol.u[0, :] * 180 / np.pi, label='No Path Constraint', color=color)

    ax_beta.plot(sol.t, sol.u[1, :] * 180 / np.pi, label='No Path Constraint', color=color)


ax_ground.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_ground.set_ylabel(r'Crossrange, $\theta$ [deg]')

ax_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
ax_hv.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

ax_gam.set_xlabel(r'Time $t$ [s]')
ax_gam.set_ylabel(r'FPA, $\gamma$ [deg]')

ax_cost.plot(gamma_0, cross_f)
ax_cost.set_xlabel(r'Initial FPA, $\gamma_0$ [deg]')
ax_cost.set_ylabel('Terminal Crossrange,\n'+r'$\theta_f$ [deg]')

ax_alpha.set_xlabel(r'Time $t$ [s]')
ax_alpha.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

ax_beta.set_xlabel(r'Time $t$ [s]')
ax_beta.set_ylabel(r'Bank Angle, $\beta$ [deg]')

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min(gamma_0), max(gamma_0)), cmap=cmap),
             cax=cax, label=r'Initial FPA, $\gamma_0$ [deg]', orientation='horizontal', aspect=50)

fig.subplots_adjust(
        top=0.905,
        bottom=0.087,
        left=0.112,
        right=0.972,
        hspace=0.659,
        wspace=0.28
)

plt.show()
