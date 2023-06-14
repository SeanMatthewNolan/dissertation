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

fig_heat = plt.figure(figsize=(3, 2.9))
ax_heat = curvelinear(fig_heat)

directory = 'star_cav_h_45_00'
files = [f'./{directory}/{file}'
         for file in os.listdir(f'./{directory}') if re.match(r'(ray)(\d*)(.bin)', file)]

data = [giuseppe.load_sol_set(file) for file in files[::]]

all_sols = []
for sols in data:
    all_sols += sols

print(len(all_sols))

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

    segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)

    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(norms)
    line = ax_heat.add_collection(lc)


cbar = fig_heat.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap), label='Terminal Velocity [ft/s]', ax=ax_heat,
        orientation='horizontal')

# ax_heat.set_title('Max Terminal Velocities\nWith Bank Angle Control')
ax_heat.set_xlabel(r'Azimuth [deg]')
ax_heat.set_ylabel(r'Range [NM]')
ax_heat.autoscale()

fig_heat.subplots_adjust(
        top=0.948,
        bottom=0.128,
        left=0.116,
        right=0.95,
        hspace=0.2,
        wspace=0.2
)

plt.show()
