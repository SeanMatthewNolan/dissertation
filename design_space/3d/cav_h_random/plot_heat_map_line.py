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
from matplotlib import rcParams, gridspec, rc

from misc import circle_ang_dist, calc_bearing, curvelinear

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


os.chdir(os.path.dirname(__file__))

SKIP = 25
# RE = 20_902_900 / 5280
RE = 20_902_900 / 6076.118

rad2deg = 180/3.141592653589793

INCLUDE_TARGETS = True
INCLUDE_FOOTPRINT = True


fig = plt.figure(figsize=(6.5, 7.5))

# gs = gridspec.GridSpec(4, 2, height_ratios=(2, 0.1, 1, 1))
#
# ax_samples = fig.add_subplot(gs[0, :], projection='3d')
# ax_path = fig.add_subplot(gs[2, 0])
# ax_hv = fig.add_subplot(gs[2, 1])
# ax_gam = fig.add_subplot(gs[3, 0])
# ax_alpha = fig.add_subplot(gs[3, 1])
# cb_ax = fig.add_subplot(gs[1, :])

gs = gridspec.GridSpec(4, 3, height_ratios=(1, 1, 1, 0.1))

ax_samples = fig.add_subplot(gs[0:2, 0:2], projection='3d')
ax_path = fig.add_subplot(gs[0, 2])
ax_hv = fig.add_subplot(gs[1, 2])
# ax_alpha = fig.add_subplot(gs[2, 0])
# ax_gam = fig.add_subplot(gs[2, 1])

ax_scatter_h = fig.add_subplot(gs[2, 0])
ax_scatter_v = fig.add_subplot(gs[2, 1])
ax_scatter_gam = fig.add_subplot(gs[2, 2])

cb_ax = fig.add_subplot(gs[3, :])

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

for sol in random.sample(sols, 250):
        color = cmap(sol.x[3, -1] / max_vf)
        ax_path.plot(np.rad2deg(sol.x[1, :]), sol.x[0, :] / 1e5, color=color)
        ax_hv.plot(sol.x[3, :] / 1e3, sol.x[0, :] / 1e5, color=color)
        # ax_gam.plot(sol.t, np.rad2deg(sol.x[4, :]), color=color)
        # ax_alpha.plot(sol.t, np.rad2deg(sol.u[0, :]), color=color)

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

    ax_scatter_h.add_collection(generate_line_collection(h0s / 1e5, vfs/1e3, vfs))
    ax_scatter_v.add_collection(generate_line_collection(v0s / 1e3, vfs/1e3, vfs))
    ax_scatter_gam.add_collection(generate_line_collection(np.rad2deg(gam0s), vfs/1e3, vfs))

    # ax_scatter_h.scatter(h0s / 1e5, vfs/1e3, c=vfs / max_vf, cmap=cmap)
    # ax_scatter_v.scatter(v0s / 1e3, vfs/1e3, c=vfs / max_vf, cmap=cmap)
    # ax_scatter_gam.scatter(np.rad2deg(gam0s), vfs/1e3, c=vfs / max_vf, cmap=cmap)

ax_samples.set_xlim(min_v0, max_v0)
ax_samples.set_ylim(min_gam0, max_gam0)
ax_samples.set_zlim(min_h0, max_h0)

ax_scatter_h.autoscale()
ax_scatter_v.autoscale()
ax_scatter_gam.autoscale()

fig.suptitle('CAV-H Maximum Terminal Velocity\nVariation of Initial States in 3D Star')

ax_samples.set_title('Sample Points')
ax_samples.set_xlabel(r'Initial Velocity, $v_0$ [$\times 10^3$ ft/s]')
ax_samples.set_ylabel(r'Initial FPA, $\gamma_0$ [deg]')
ax_samples.set_zlabel('Initial Altitude,' + r'$h_0$ [$\times 10^5$ ft]')

ax_path.set_title('Flight Path')
ax_path.set_xlabel(r'Downrange, $\phi_f$ [deg]')
ax_path.set_ylabel(r'Altitude, $h$ [$\times 10^5$ ft]')

ax_hv.set_title('h-v Plot')
ax_hv.set_xlabel(r'Velocity, $v$ [$\times 10^3$ ft/s]')
ax_hv.set_ylabel(r'Altitude, $h$ [$\times 10^5$ ft]')

ax_scatter_h.set_title('Initial Altitude Effect')
ax_scatter_h.set_xlabel('Initial Altitude,\n' + r'$h_0$ [$\times 10^5$ ft]')
ax_scatter_h.set_ylabel('Terminal Velocity,\n' + r'$v_f$ [$\times 10^3$ ft/s]')

ax_scatter_v.set_title('Initial Velocity Effect')
ax_scatter_v.set_xlabel('Initial Velocity,\n' + r'$v_0$ [$\times 10^3$ ft/s]')
ax_scatter_v.set_ylabel('Terminal Velocity,\n' + r'$v_f$ [$\times 10^3$ ft/s]')

ax_scatter_gam.set_title('Initial FPA Effect')
ax_scatter_gam.set_xlabel('Initial FPA,\n' + r'$\gamma_0$ [deg]')
ax_scatter_gam.set_ylabel('Terminal Velocity,\n' + r'$v_f$ [$\times 10^3$ ft/s]')

# ax_alpha.set_title('Angle-of-Attack')
# ax_alpha.set_xlabel(r'Time, $t$, [sec]')
# ax_alpha.set_ylabel(r'AoA, $\alpha$ [deg]')
#
# ax_gam.set_title('Flight Path Angle')
# ax_gam.set_xlabel(r'Time, $t$, [sec]')
# ax_gam.set_ylabel(r'FPA, $\gamma$ [deg]')

cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        orientation='horizontal', label='Terminal Velocity [ft/s]', cax=cb_ax)

fig.subplots_adjust(
        top=0.883,
        bottom=0.078,
        left=0.129,
        right=0.943,
        hspace=0.963,
        wspace=0.774
)

plt.show()
