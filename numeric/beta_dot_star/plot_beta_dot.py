import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection

from misc import curvelinear, circle_ang_dist, calc_bearing

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm,siunitx}'

RE = 20_902_900 / 6076.118

sol_set = gp.load_sol_set('beta_dot_ref.data')
footprint_set = gp.load_sol_set('cav_h_footprint.bin')
phi_f0, theta_f0 = sol_set[0].x[1, -1], sol_set[0].x[2, -1]

norm_vals = np.array(
        [circle_ang_dist(phi_f0, theta_f0, sol.x[1, -1], sol.x[2, -1])
         for sol in sol_set]) * RE

cmap = plt.colormaps['viridis']
norm = Normalize(min(norm_vals), max(norm_vals))

fig = plt.figure(figsize=(6.5, 5))
title = fig.suptitle(r'HGV Maximum Terminal Velocity Trajectories: Bank Rate Control')

gs = gridspec.GridSpec(3, 2, height_ratios=[1, 1, 0.1])

ax_foot = curvelinear(fig, gs=gs[1, 0])
ax_lift = fig.add_subplot(gs[1, 1])

ax_beta = fig.add_subplot(gs[0, 0])
ax_beta_dot = fig.add_subplot(gs[0, 1])

ax_cb = fig.add_subplot(gs[-1, :])

c_l0, c_l1, c_l2, c_l3 = sol_set[-1].k[4:8]
a = sol_set[-1].k[12]


def compute_cl(alpha, v):
    return c_l1 * alpha + c_l2 * np.exp(c_l3 * v / a) + c_l0


phi_fs, theta_fs = [], []
for sol, norm_val in zip(sol_set, norm_vals):
    color = cmap(norm(norm_val))

    eigs = np.array([
        np.linalg.eig(_h_uu)[0] for _h_uu in sol.aux['H_uu']
    ])

    # ax_foot.plot(sol.x[1, :] * 180 / np.pi, sol.x[2, :] * 180 / np.pi, sol.x[0, :] / 1e5, color=color)
    ax_lift.plot(sol.t, compute_cl(sol.u[0, :], sol.x[3, :]), color=color)
    ax_beta.plot(sol.t, np.rad2deg((sol.x[6, :] + np.pi) % (2 * np.pi) - np.pi), color=color)
    # ax_alpha.plot(sol.t, eigs[:, 1], color=color)
    ax_beta_dot.plot(sol.t, sol.u[1, :] * 180 / np.pi, color=color)

    phi_fs.append(sol.x[1, -1])
    theta_fs.append(sol.x[2, -1])

star_theta = np.array(theta_fs)
star_phi = np.array(phi_fs)

r = RE * circle_ang_dist(0., 0., star_phi, star_theta)
bear = calc_bearing(0., 0., star_phi, star_theta)

x = r * np.cos(bear)
y = r * np.sin(bear)

points = np.array([x, y]).T.reshape(-1, 1, 2)

segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)

footprint_phi = []
footprint_theta = []
for sol in footprint_set:
    footprint_phi.append(sol.x[1, -1])
    footprint_theta.append(sol.x[2, -1])

footprint_phi = np.array(footprint_phi)
footprint_theta = np.array(footprint_theta)

footprint_range = RE * circle_ang_dist(0., 0., footprint_phi, footprint_theta)
footprint_bear = calc_bearing(0., 0., footprint_phi, footprint_theta)

ax_foot.plot(footprint_range * np.cos(footprint_bear), footprint_range * np.sin(footprint_bear), ':',
             color='C7', label='Reachable Footprint')

ax_foot.legend(loc='upper center', fontsize=6)

lc = LineCollection(segments, cmap=cmap, norm=norm)
lc.set_array(norm_vals)
line = ax_foot.add_collection(lc)

fig.colorbar(plt.cm.ScalarMappable(norm=norm), cmap=cmap,
             cax=ax_cb, label=r'Distance Traveled in Continuation [NM]',
             orientation='horizontal')

ax_beta.set_xlabel(r'Time, $t$ [s]')
ax_beta.set_ylabel(r'Bank, $\beta$ [deg]')

ax_beta_dot.set_xlabel(r'Time, $t$ [s]')
ax_beta_dot.set_ylabel(r'Bank Rate, $\dot{\beta}$ [deg/s]')

ax_lift.set_xlabel(r'Time, $t$ [s]')
ax_lift.set_ylabel(r'Lift Coeff., $C_l$')

ax_foot.set_xlabel(r'Azimuth [deg]')
ax_foot.set_ylabel(r'Range [NM]')
ax_foot.autoscale()

fig.subplots_adjust(
        top=0.91,
        bottom=0.113,
        left=0.116,
        right=0.97,
        hspace=0.592,
        wspace=0.285
)

plt.show()
