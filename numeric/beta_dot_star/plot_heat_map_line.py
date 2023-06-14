import os.path
import re
import random

import numpy as np
from matplotlib import rcParams, rc, gridspec
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LogNorm
from scipy.integrate import simpson

import giuseppe

from misc import circle_ang_dist, calc_bearing, curvelinear


os.chdir(os.path.dirname(__file__))

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

SKIP = 1
RE = 20_902_900 / 6076.118

rad2deg = 180/3.141592653589793

vf_range = (5000, 15000)
# vf_range = (1e-6, 1)

cmap = plt.colormaps['viridis']
norm = Normalize(*vf_range)
# norm = LogNorm(*vf_range)

fig_traj = plt.figure(figsize=(4.5, 4.5))
fig_heat = plt.figure(figsize=(6.5, 4.5))
ax_traj = fig_traj.add_subplot(projection='3d')
ax_heat = curvelinear(fig_heat)

# gs = gridspec.GridSpec(2, 1)
#
# fig = plt.figure(figsize=(6.5, 7.5))
# ax_traj = fig.add_subplot(gs[0, 0], projection='3d')
# ax_heat = curvelinear(fig, gs[1, 0])

directory = 'star_cav_h_45_00'
files = [f'./{directory}/{file}'
         for file in os.listdir(f'./{directory}') if re.match(r'(ray)(\d*)(.bin)', file)]

data = [giuseppe.load_sol_set(file) for file in files[::]]

all_sols = []
for sols in data:
    all_sols += sols

print(len(all_sols))

for sol in random.sample(all_sols, 100):
    ax_traj.plot(np.rad2deg(sol.x[1, :]), np.rad2deg(sol.x[2, :]), sol.x[0, :] / 1000, color='C0', linewidth=0.5)


for sols in data:
    star_theta = []
    star_phi = []
    norms = []

    for sol in sols[::SKIP]:
        theta = sol.x[2, :]
        phi = sol.x[1, :]
        # norms.append(sol.x[3, -1] - 5_000) / (15_000 - 5_000)
        norms.append(sol.x[3, -1])
        # norms.append(np.min(sol.aux['H_uu'][:, 1, 1]))

        star_theta.append(theta[-1])
        star_phi.append(phi[-1])

    star_theta = np.array(star_theta)
    star_phi = np.array(star_phi)

    r = RE * circle_ang_dist(0., 0., star_phi, star_theta)
    bear = calc_bearing(0., 0., star_phi, star_theta)

    x = r * np.cos(bear)
    y = r * np.sin(bear)

    points = np.array([x, y]).T.reshape(-1, 1, 2)

    ax_traj.plot(np.rad2deg(star_phi), np.rad2deg(star_theta), 0, color='C1', zorder=0, linewidth=0.5)

    segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)

    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(norms)
    line = ax_heat.add_collection(lc)


cbar = fig_heat.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap), label='Terminal Velocity [ft/s]', ax=ax_heat)

# sols_footprint = giuseppe.load_sol_set('cav_h_footprint.bin')
#
# footprint_phi = []
# footprint_theta = []
# for sol in sols_footprint:
#     footprint_phi.append(sol.x[1, -1])
#     footprint_theta.append(sol.x[2, -1])
#
# footprint_phi = np.array(footprint_phi)
# footprint_theta = np.array(footprint_theta)
#
# footprint_range = RE * circle_ang_dist(0., 0., footprint_phi, footprint_theta)
# footprint_bear = calc_bearing(0., 0., footprint_phi, footprint_theta)
#
# ax_heat.plot(footprint_range * np.cos(footprint_bear), footprint_range * np.sin(footprint_bear), ':',
#         color='C7', label='Reachable Footprint')

# ax_heat.legend(loc='upper center')

ax_traj.set_title('Maximum Terminal Velocity Trajectories\nWith Bank Angle Control')
ax_traj.set_xlabel(r'Downrange [deg]')
ax_traj.set_ylabel(r'Crossrange [deg]')
ax_traj.set_zlabel(r'Altitude [$\times 10^5$ ft]')

ax_heat.set_title('Maximum Terminal Velocity Trajectories\nWith Bank Angle Control')
ax_heat.set_xlabel(r'Azimuth [deg]')
ax_heat.set_ylabel(r'Range [NM]')
ax_heat.autoscale()

plt.show()
