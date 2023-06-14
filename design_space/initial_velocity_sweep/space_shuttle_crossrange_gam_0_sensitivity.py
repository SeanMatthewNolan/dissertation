import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize


rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

fpa_sweep = gp.load_sol_set('initial_fpa_sweep.data')

fig = plt.figure(figsize=(6.5, 5))
title = fig.suptitle('Space Shuttle Maximum Crossrange\n'
                     'Demonstration of Adjoint Variable Relationship to Cost Functional')

gs = gridspec.GridSpec(3, 2, height_ratios=(1, 0.05, 1))

ax_ground = fig.add_subplot(gs[0, 1])
ax_hv = fig.add_subplot(gs[0, 0])

cax = fig.add_subplot(gs[1, :])

ax_cost = fig.add_subplot(gs[2, 0])
ax_diff = fig.add_subplot(gs[2, 1])


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
nu_gamma = []
for sol in fpa_sweep:
    gamma_0.append(np.rad2deg(sol.x[4, 0]))
    cross_f.append(np.rad2deg(sol.x[2, -1]))
    nu_gamma.append(sol.nu0[5])

for sol in fpa_sweep[::-1]:
    color = cmap((np.rad2deg(sol.x[4, 0]) - min(gamma_0)) / (max(gamma_0) - min(gamma_0)))

    ax_ground.plot(sol.x[1, :] * 180 / np.pi, sol.x[2, :] * 180 / np.pi, color=color)
    ax_hv.plot(sol.x[3, :] / 1e3, sol.x[0, :] / 1e5, color=color)


ax_ground.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_ground.set_ylabel(r'Crossrange, $\theta$ [deg]')

ax_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
ax_hv.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

ax_cost.plot(gamma_0, cross_f, linewidth=2)
ax_cost.set_xlabel(r'Initial FPA, $\gamma_0$ [deg]')
ax_cost.set_ylabel('Terminal Crossrange,\n'+r'$\theta_f$ [deg]')

ax_diff.plot((np.asarray(gamma_0[:-1]) + np.asarray(gamma_0[1:]))/2, np.diff(cross_f) / np.diff(gamma_0),
             linewidth=2, label='From Continuation')
ax_diff.plot(gamma_0, nu_gamma, ':', linewidth=2, label=r'$-\nu_{\gamma_0}$')
ax_diff.set_xlabel(r'Initial FPA, $\gamma_0$ [deg]')
ax_diff.set_ylabel(r'Cost Sensitivity $\frac{d \theta_f}{d \gamma_0}$')
ax_diff.legend()

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min(gamma_0), max(gamma_0)), cmap=cmap),
             cax=cax, label=r'Initial FPA, $\gamma_0$ [deg]', orientation='horizontal', aspect=50)

fig.subplots_adjust(
        top=0.873,
        bottom=0.113,
        left=0.112,
        right=0.972,
        hspace=0.664,
        wspace=0.279
)

plt.show()
