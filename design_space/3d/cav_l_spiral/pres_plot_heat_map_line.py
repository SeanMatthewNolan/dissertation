import os
import os.path
import re
import pickle
import random

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.colors import Normalize
from matplotlib import rc, rcParams, gridspec

from misc import circle_ang_dist, calc_bearing, curvelinear

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


os.chdir(os.path.dirname(__file__))

SKIP = 25
# RE = 20_902_900 / 5280
RE = 20_902_900 / 6076.118

rad2deg = 180/3.141592653589793

fig = plt.figure(figsize=(3, 2.9))

ax_samples = fig.add_subplot(projection='3d')

max_vf = 24_000

cmap = plt.colormaps['viridis']
norm = Normalize(0, max_vf)


def generate_line_collection(_x, _y, _norms):
    points = np.array([_x, _y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(_norms)
    return lc


def generate_line_collection_3d(_x, _y, _z, _norms):
    points = np.array([_x, _y, _z]).T.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1, :], points[1:, :]], axis=1)
    lc = Line3DCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(_norms)
    return lc


directory = 'data'
sols = []
files = []
files += [f'./{directory}/{file}'
          for file in os.listdir(f'./{directory}') if re.match(r'(ray)(\d*)(.bin)', file)]

sol_series = []
for file in files[::]:
    # print(file)
    with open(file, 'rb') as f:
        sols_i = pickle.load(f)
        sols += sols_i.solutions
        sol_series.append(sols_i.solutions)

min_h0 = 1e99
min_v0 = 1e99
min_gam0 = 1e99

max_h0 = 0
max_v0 = 0
max_gam0 = 0

sample_lines = []
for sols in random.sample(sol_series, len(sol_series)):
    h0_list = []
    v0_list = []
    gam0_list = []
    vf_list = []
    for sol in sols:
        h0_list.append(sol.x[0, 0])
        v0_list.append(sol.x[3, 0])
        gam0_list.append(sol.x[4, 0])
        vf_list.append(sol.x[3, -1])

    min_h0 = min(min(h0_list)/1e5, min_h0)
    min_v0 = min(min(v0_list)/1e3, min_v0)
    min_gam0 = min(np.rad2deg(min(gam0_list)), min_gam0)

    max_h0 = max(max(h0_list)/1e5, max_h0)
    max_v0 = max(max(v0_list)/1e3, max_v0)
    max_gam0 = max(np.rad2deg(max(gam0_list)), max_gam0)

    h0s = np.array(h0_list)
    v0s = np.array(v0_list)
    gam0s = np.array(gam0_list)
    vfs = np.array(vf_list)

    # ax_samples.scatter(v0s / 1e3, np.rad2deg(gam0s), h0s / 1e5, c=vfs / max_vf, cmap=cmap)
    sample_lines.append(
            ax_samples.add_collection(generate_line_collection_3d(v0s / 1e3, np.rad2deg(gam0s), h0s / 1e5, vfs)))


ax_samples.set_xlim(min_v0, max_v0)
ax_samples.set_ylim(min_gam0, max_gam0)
ax_samples.set_zlim(min_h0, max_h0)

ax_samples.set_title(r'CAV-L Max $v_f$:' '\nVariation of Initial States')
ax_samples.set_xlabel(r'$v_0$ [$\times 10^3$ ft/s]')
ax_samples.set_ylabel(r'$\gamma_0$ [deg]')
ax_samples.set_zlabel(r'$h_0$ [$\times 10^5$ ft]')

cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        orientation='horizontal', label='Terminal Velocity [ft/s]', ax=ax_samples, aspect=50)

fig.subplots_adjust(
        top=0.817,
        bottom=0.095,
        left=0.062,
        right=0.95,
        hspace=0.963,
        wspace=0.774
)

plt.show()
