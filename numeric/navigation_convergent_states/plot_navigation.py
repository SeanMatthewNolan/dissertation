import pickle

import matplotlib.pyplot as plt
import giuseppe as gp
import numpy as np
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import LogNorm


rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

cmap = plt.colormaps['viridis']
norm = LogNorm(0.01, 0.1)

fig_ex = plt.figure(figsize=(6.5, 7.5))
fig_ex.suptitle('Navigation Cart Problem: Variation in ' + r'$\sigma_\omega$')

gs = gridspec.GridSpec(5, 2, height_ratios=[1, 1, 1, 1, 0.1])

ax_path_path = fig_ex.add_subplot(gs[0, 0])
# ax_path_path.plot(500, 250, '*', label='Beacon', color='C3')
ax_path_path.set_title('Position: Path Cost')
ax_path_path.set_xlabel(r'Position, $x$ [ft]')
ax_path_path.set_ylabel(r'Position, $y$ [ft]')

ax_theta_path = fig_ex.add_subplot(gs[1, 0])
ax_theta_path.set_title('Heading: Path Cost')
ax_theta_path.set_xlabel(r'$t$ [s]')
ax_theta_path.set_ylabel(r'Heading Angle,' + '\n' + r'$\theta$ [deg]')

ax_omega_path = fig_ex.add_subplot(gs[2, 0])
ax_omega_path.set_title('Turn-Rate: Path Cost')
ax_omega_path.set_xlabel(r'$t$ [s]')
ax_omega_path.set_ylabel(r'Turn Rate,' + '\n' + r'$\omega$ [deg / s]')

ax_cost_path = fig_ex.add_subplot(gs[3, 0])
ax_cost_path.set_title('Variance: Path Cost')
ax_cost_path.set_xlabel(r'$t$ [s]')
ax_cost_path.set_ylabel('Position Variance,\n' + r'$p_{xx} + p_{yy}$ [$\text{ft}^2$]')

with open('path.data', 'rb') as file:
    path_data = pickle.load(file)

for sol in path_data['final_sols'][::-1]:
    xb_path, yb_path = sol.k[1], sol.k[2]
    sig_omega = sol.k[4]
    color = cmap(norm(sig_omega))

    ax_path_path.plot(sol.x[0, :], sol.x[1, :], color=color)
    ax_theta_path.plot(sol.t, np.rad2deg(sol.x[2, :]), color=color)
    ax_omega_path.plot(sol.t, np.rad2deg(sol.u[0, :]), color=color)
    ax_cost_path.plot(sol.t, sol.x[3, :] + sol.x[6, :], color=color)


ax_path_term = fig_ex.add_subplot(gs[0, 1])
# ax_path_term.plot(500, 250, '*', label='Beacon', color='C3')
ax_path_term.set_title('Position: Terminal Cost')
ax_path_term.set_xlabel(r'Position, $x$ [ft]')
ax_path_term.set_ylabel(r'Position, $y$ [ft]')

ax_theta_term = fig_ex.add_subplot(gs[1, 1])
ax_theta_term.set_title('Heading: Terminal Cost')
ax_theta_term.set_xlabel(r'$t$ [s]')
ax_theta_term.set_ylabel(r'Heading Angle,' + '\n' + r'$\theta$ [deg]')

ax_omega_term = fig_ex.add_subplot(gs[2, 1])
ax_omega_term.set_title('Turn-Rate: Terminal Cost')
ax_omega_term.set_xlabel(r'$t$ [s]')
ax_omega_term.set_ylabel(r'Turn Rate,' + '\n' + r'$\omega$ [deg / s]')

ax_cost_term = fig_ex.add_subplot(gs[3, 1])
ax_cost_term.set_title('Variance: Terminal Cost')
ax_cost_term.set_xlabel(r'$t$ [s]')
ax_cost_term.set_ylabel('Position Variance,\n' + r'$p_{xx} + p_{yy}$ [$\text{ft}^2$]')

with open('term.data', 'rb') as file:
    term_data = pickle.load(file)

for sol in term_data['final_sols'][::-1]:
    xb_term, yb_term = sol.k[1], sol.k[2]
    sig_omega = sol.k[4]
    color = cmap(norm(sig_omega))

    ax_path_term.plot(sol.x[0, :], sol.x[1, :], color=color)
    ax_theta_term.plot(sol.t, np.rad2deg(sol.x[2, :]), color=color)
    ax_omega_term.plot(sol.t, np.rad2deg(sol.u[0, :]), color=color)
    ax_cost_term.plot(sol.t, sol.x[3, :] + sol.x[6, :], color=color)

ax_cb = fig_ex.add_subplot(gs[-1, :])
fig_ex.colorbar(plt.cm.ScalarMappable(norm=norm), cmap=cmap,
                cax=ax_cb, label=r'Turn-Rate Standard Deviation $\sigma_\omega$, [deg/s]',
                orientation='horizontal', aspect=50)
fig_ex.tight_layout()

fig_ex.subplots_adjust(
        top=0.914,
        bottom=0.075,
        left=0.128,
        right=0.955,
        hspace=1.0,
        wspace=0.365
)

plt.show()
