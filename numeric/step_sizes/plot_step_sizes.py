import re
import os
import csv

import giuseppe
import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec,rc
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
cmap = plt.colormaps['viridis']

fig = plt.figure(figsize=(6.5, 7.5))
title = fig.suptitle('Effect of Number of Steps on Bisection Method')

gs = gridspec.GridSpec(4, 2, height_ratios=(1, 1, 1, 0.1))

ax_time = fig.add_subplot(gs[0, 0])
ax_num_bisection = fig.add_subplot(gs[0, 1])

# ax_traj_050 = fig.add_subplot(gs[1, 0])
# ax_traj_100 = fig.add_subplot(gs[1, 1])
# ax_traj_250 = fig.add_subplot(gs[2, 0])
# ax_traj_500 = fig.add_subplot(gs[2, 1])

ax_cb = fig.add_subplot(gs[3, :])

norm = Normalize(10_000, 17_500)
vs = []

for idx, num_steps in enumerate([25, 50, 100, 250]):
    ax = fig.add_subplot(gs[idx + 2])

    file_name = f'cav_h_{num_steps:d}.bin'
    sol_set = gp.load_sol_set(file_name)

    for sol in sol_set:
        ax.plot(np.rad2deg(sol.x[1, :]), sol.x[0, :] / 1e5, color=cmap(norm(sol.x[2, -1])), linewidth=0.5)
        vs.append(sol.x[2, -1])

    ax.set_title(f'User Specified {num_steps:d} Steps')
    ax.set_xlabel('Downrange Angle, $\phi$ [deg]')
    ax.set_ylabel('Altitude,\n' + r'$h$ [$\times 10^5$ ft]')


fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_cb,
             orientation='horizontal', label='Constrained Terminal Downrange [deg]')


with open('times.csv', 'r') as csv_file:
    rows = [row for row in csv.reader(csv_file)]
    num_steps = np.round(np.asarray(rows[0], dtype=float))
    times = np.asarray(rows[1], dtype=float) * 1e-9
    num_bisects = np.round(np.asarray(rows[2], dtype=float))

# ax_time.plot(num_steps, times, '*--')
ax_time.plot(num_steps, times)
# ax_time.set_xscale('log')
# ax_time.set_yscale('log')
ax_time.set_title(f'Continuation Times')
ax_time.set_xlabel('Number of User Specified Steps')
ax_time.set_ylabel('Continuation Time [s]')

# ax_num_bisection.plot(num_steps, num_bisects, '+--', markersize=2)
ax_num_bisection.plot(num_steps, num_bisects)
# ax_num_bisection.set_xscale('log')
# ax_num_bisection.set_yscale('log')
ax_num_bisection.set_title(f'Number of Bisections')
ax_num_bisection.set_xlabel('Number of User Specified Steps')
ax_num_bisection.set_ylabel('Number of Bisections')

print(min(vs))
print(max(vs))

fig.subplots_adjust(
        top=0.914,
        bottom=0.075,
        left=0.102,
        right=0.977,
        hspace=0.763,
        wspace=0.265
)

plt.show()
