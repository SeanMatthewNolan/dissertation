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
          for file in os.listdir(f'./') if re.match(r'(vf_first_)(\d*)(.bin)', file)]

data = [giuseppe.load_sol_set(file) for file in files[::]]

cmap = plt.colormaps['viridis']

for sols in data:
    color = cmap((sols[-1].x[3, -1] - 5_000) / (15_000 - 5_000))

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

# ax.plot(0, 0, label='Initial Point', marker='+', linestyle='')

ax.set_title('Reachability Footprints\nas Function of Terminal Velocity')
ax.set_xlabel(r'Angle to Glide Origin')
ax.set_ylabel(r'Range to Glide Origin [NM]')

ax.autoscale()

# ax.legend(title='Time [s]', loc='upper left', bbox_to_anchor=(1, 1), fontsize=10, title_fontsize=10)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(5_000, 15_000), cmap=cmap), ax=ax, label=r'Terminal Velocity [ft/s]')
# fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(5_000, 15_000), cmap=cmap), ax=ax, orientation='horizontal',
#              aspect=50, label=r'Terminal Velocity [ft/s]')

fig.subplots_adjust(
        top=0.897,
        bottom=0.096,
        left=0.105,
        right=0.977,
        hspace=0.2,
        wspace=0.2
)

plt.show()
