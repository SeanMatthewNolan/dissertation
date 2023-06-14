import os
import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec

import giuseppe.data_classes.solution_set
from misc import circle_ang_dist, calc_bearing, curvelinear, load

os.chdir(os.path.dirname(__file__))

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

SKIP = 1
# RE = 20_902_900 / 5280
RE = 20_902_900 / 6076.118
rad2deg = 180 / 3.141592653589793

left_sols = giuseppe.data_classes.load_sol_set('cav_l_footprint_left.bin')
right_sols = giuseppe.data_classes.load_sol_set('cav_l_footprint_right.bin')

fig = plt.figure(figsize=(6.5, 5))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1])

ax = curvelinear(fig, gs[1, 0])

ax_traj = fig.add_subplot(gs[0, 0], projection='3d')
ax_alpha = fig.add_subplot(gs[0, 1])
ax_bank = fig.add_subplot(gs[1, 1])

xi_select = np.deg2rad(230)

vfss = []
for sols, label, color in ((left_sols, 'Left', 'C0'), (right_sols, 'Right', 'C1')):
    vfs = []
    footprint_phi = []
    footprint_theta = []

    for sol in sols[::SKIP]:
        phi = sol.x[1, :]
        theta = sol.x[2, :]
        xi = sol.k[21]

        if np.isclose(xi % (2 * np.pi), xi_select):
            r_i = RE * circle_ang_dist(0., 0., phi, theta)
            bear_i = calc_bearing(0., 0., phi, theta)
            ax.plot(r_i * np.cos(bear_i), r_i * np.sin(bear_i), color=color,
                    label=label + r': $\xi$ = ' + f'{np.rad2deg(xi):3.0f} deg')
            ax_traj.plot(np.rad2deg(theta), np.rad2deg(phi), sol.x[0, :] / 10_000, color=color)
            ax_alpha.plot(sol.t, np.rad2deg(sol.u[0, :]), color=color)
            ax_bank.plot(sol.t, np.rad2deg(sol.u[1, :]), color=color)

        footprint_phi.append(phi[-1])
        footprint_theta.append(theta[-1])
        vfs.append(sol.x[3, -1])

    footprint_phi = np.array(footprint_phi)
    footprint_theta = np.array(footprint_theta)

    r = RE * circle_ang_dist(0., 0., footprint_phi, footprint_theta)
    bear = calc_bearing(0., 0., footprint_phi, footprint_theta)

    ax.plot(r * np.cos(bear), r * np.sin(bear), '--',  color=color, label=label + ': Footprint')
    vfss.append(vfs)


fig.suptitle('Reachability Footprint: Multiple Solutions')

ax.set_xlabel(r'Angle to Glide Origin [deg]')
ax.set_ylabel(r'Range to Glide Origin [NM]')

ax_traj.set_xlabel(r'Crossrange [deg]')
ax_traj.set_ylabel(r'Downrange [deg]')
ax_traj.set_zlabel(r'Altitude [$\times 10^5$ ft]')

ax_alpha.set_xlabel(r'Time [sec]')
ax_alpha.set_ylabel(r'AOA [deg]')

ax_bank.set_xlabel(r'Time [sec]')
ax_bank.set_ylabel(r'Bank [deg]')

fig.legend(ncols=2, loc='lower center')

ax.autoscale()

fig.subplots_adjust(
        top=0.876,
        bottom=0.203,
        left=0.055,
        right=0.964,
        hspace=0.477,
        wspace=0.157
)

plt.show()
