import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

fpa_sweep = gp.load_sol_set('initial_fpa_sweep.data')

fig = plt.figure(figsize=(3.5, 2.9))
title = fig.suptitle('Effect of Initial FPA on Crossrange')

gs = gridspec.GridSpec(2, 1)

ax_cost = fig.add_subplot(gs[1, 0])
ax_ground = fig.add_subplot(gs[0, 0])


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
    ax_ground.plot(sol.x[1, :] * 180 / np.pi, sol.x[2, :] * 180 / np.pi, color=color)


ax_ground.set_xlabel(r'$\phi$ [deg]')
ax_ground.set_ylabel(r'$\theta$ [deg]')

ax_cost.plot(gamma_0, cross_f)
ax_cost.set_xlabel(r'$\gamma_0$ [deg]')
ax_cost.set_ylabel(r'$\theta_f$ [deg]')

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min(gamma_0), max(gamma_0)), cmap=cmap),
             ax=ax_ground, label=r'$\gamma_0$ [deg]', orientation='vertical')

fig.subplots_adjust(
        top=0.884,
        bottom=0.195,
        left=0.198,
        right=0.9,
        hspace=0.7,
        wspace=0.28
)

plt.show()
