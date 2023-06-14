import os
import os.path
import re
import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams, rc
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

import giuseppe

from misc import circle_ang_dist, calc_bearing, curvelinear


os.chdir(os.path.dirname(__file__))

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

SKIP = 1
# RE = 20_902_900 / 5280
RE = 20_902_900 / 6076.118

rad2deg = 180/3.141592653589793

INCLUDE_TARGETS = True
INCLUDE_FOOTPRINT = True

vf_range = (1000, 15_000)

cmap = plt.colormaps['viridis']
norm = Normalize(*vf_range)

fig = plt.figure(figsize=(6.5, 4.5))
ax = curvelinear(fig)

sols = []
files = []
for directory in tuple(os.walk('.'))[0][1]:
    if re.match(r'(star_cav_h_)(\d\d)(_)(-?)(\d\d)', directory):
        # print(directory)
        files += [f'./{directory}/{file}'
                  for file in os.listdir(f'./{directory}') if re.match(r'(ray)(\d*)(.bin)', file)]


data = [giuseppe.load_sol_set(file) for file in files[::]]

for sols in data:
    footprint_theta = []
    footprint_phi = []
    norms = []

    for sol in sols[::SKIP]:
        theta = sol.x[1, :]
        phi = sol.x[2, :]
        # norms.append(sol.x[3, -1] - 5_000) / (15_000 - 5_000)
        norms.append(sol.x[3, -1])

        footprint_theta.append(theta[-1])
        footprint_phi.append(phi[-1])

    footprint_theta = np.array(footprint_theta)
    footprint_phi = np.array(footprint_phi)

    r = RE * circle_ang_dist(0., 0., footprint_theta, footprint_phi)
    bear = calc_bearing(0., 0., footprint_theta, footprint_phi)

    x = r * np.cos(bear)
    y = r * np.sin(bear)

    points = np.array([x, y]).T.reshape(-1, 1, 2)

    segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)

    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(norms)
    line = ax.add_collection(lc)

# for _phi, _theta, v_f in zip(np.rad2deg(target_phi), np.rad2deg(target_theta), target_vf):
#     ax2.plot(_phi, _theta, marker='o', markersize=1, color=cmap(v_f/max_vf))

cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap), label='Terminal Velocity [ft/s]', ax=ax)

if INCLUDE_FOOTPRINT:
    with open('./cav_h_footprint.bin', 'rb') as f:
        sols_footprint = pickle.load(f)

    footprint_phi = []
    footprint_theta = []
    for sol in sols_footprint:
        footprint_phi.append(sol.x[1, -1])
        footprint_theta.append(sol.x[2, -1])

    footprint_phi = np.array(footprint_phi)
    footprint_theta = np.array(footprint_theta)

    footprint_range = RE * circle_ang_dist(0., 0., footprint_phi, footprint_theta)
    footprint_bear = calc_bearing(0., 0., footprint_phi, footprint_theta)

    ax.plot(footprint_range * np.cos(footprint_bear), footprint_range * np.sin(footprint_bear),
             color='C1', label='Reachable Footprint')

    ax.legend(loc='upper right', numpoints=1, fontsize=10, markerscale=10, framealpha=1)

print(min(norms))
print(max(norms))

# if INCLUDE_TARGETS:
#     ax1.legend(loc="best", numpoints=1, fontsize=10, markerscale=10, framealpha=1)

ax.set_title('Maximum Terminal Velocity Heat Map')
ax.set_xlabel(r'Azimuth [deg]')
ax.set_ylabel(r'Range [NM]')

plt.show()
