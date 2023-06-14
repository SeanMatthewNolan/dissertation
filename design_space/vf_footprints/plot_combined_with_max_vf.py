import os
import re

import giuseppe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams, rc
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

from misc import circle_ang_dist, calc_bearing, curvelinear, load

os.chdir(os.path.dirname(__file__))

SKIP = 1
RE = 20_902_900 / 6076.118

rad2deg = 180/3.141592653589793

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

fig = plt.figure(figsize=(6.5, 4.5))
ax = curvelinear(fig)

sols = []
files = []
files += [f'./{file}'
          for file in os.listdir(f'./') if re.match(r'(xi_first_)(-?)(\d*)(.bin)', file)]

data = [giuseppe.load_sol_set(file) for file in files[::]]

cmap = plt.colormaps['viridis']
norm = Normalize(000, 20_000)

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

    # ax.plot(r * np.cos(bear), r * np.sin(bear), color=np.array(norms))


files = []
files += [f'./{file}'
          for file in os.listdir(f'./') if re.match(r'(vf_first_)(\d*)(.bin)', file)]

data = [giuseppe.load_sol_set(file) for file in files[::]]

cmap = plt.colormaps['viridis']

for sols in data:
    color = cmap((sols[-1].x[3, -1] - 0) / (20_000 - 000))

    footprint_theta = []
    footprint_phi = []

    for sol in sols[::SKIP]:
        theta = sol.x[1, :]
        phi = sol.x[2, :]

        footprint_theta.append(theta[-1])
        footprint_phi.append(phi[-1])

    footprint_theta = np.array(footprint_theta)
    footprint_phi = np.array(footprint_phi)

    r = RE * circle_ang_dist(0., 0., footprint_theta, footprint_phi)
    bear = calc_bearing(0., 0., footprint_theta, footprint_phi)

    ax.plot(r * np.cos(bear), r * np.sin(bear), color=color)


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

# ax.plot(0, 0, label='Initial Point', marker='+', linestyle='')

ax.set_title('Comparison of Maximum Range and\nMaximum Terminal Energy Results')
ax.set_xlabel(r'Angle to Glide Origin')
ax.set_ylabel(r'Range to Glide Origin [NM]')

ax.autoscale()

# ax.legend(title='Time [s]', loc='upper left', bbox_to_anchor=(1, 1), fontsize=10, title_fontsize=10)
fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label=r'Terminal Velocity [ft/s]')

fig.subplots_adjust(
        top=0.896,
        bottom=0.096,
        left=0.105,
        right=0.977,
        hspace=0.2,
        wspace=0.2
)

plt.show()
