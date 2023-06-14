import re
import os
import csv

import giuseppe
import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc
from matplotlib.colors import Normalize, LogNorm

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
cmap = plt.colormaps['viridis']

fig = plt.figure(figsize=(6.5, 4))
title = fig.suptitle('Effect of Node Buffer Size on Bisection Method')

gs = gridspec.GridSpec(2, 2)

ax_time = fig.add_subplot(gs[0, 0])
ax_num_bisection = fig.add_subplot(gs[0, 1])
ax_final_nodes = fig.add_subplot(gs[1, 0])
ax_nodes = fig.add_subplot(gs[1, 1])

# ax_traj_050 = fig.add_subplot(gs[1, 0])
# ax_traj_100 = fig.add_subplot(gs[1, 1])
# ax_traj_250 = fig.add_subplot(gs[2, 0])
# ax_traj_500 = fig.add_subplot(gs[2, 1])

# ax_cb = fig.add_subplot(gs[1, :])

# norm = Normalize(5, 1000)
norm = LogNorm(5, 250)

final_nodes = []
buffer_sizes = []
max_phi_fs = np.deg2rad(60)
for file in [file for file in os.listdir('.') if re.match(r'(cav_h_buffer_)(\d+)(.bin)', file)]:
    sol_set = gp.load_sol_set(file)

    num_nodes = []
    phi_fs = []
    for sol in sol_set:
        # if sol.x[1, -1] > max_phi_fs + 1e-3:
        #     break

        num_nodes.append(len(sol.t))
        phi_fs.append(sol.x[1, -1])

    buffer_sizes.append(int(re.findall(r'(\d+)', file)[0]))

    ax_nodes.plot(np.rad2deg(phi_fs), num_nodes, color=cmap(norm(int(re.findall(r'(\d+)', file)[0]))))
    final_nodes.append(num_nodes[-1])


# fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax_cb,
#              orientation='horizontal', label='Constrained Terminal Downrange [deg]')


with open('buffers.csv', 'r') as csv_file:
    rows = [row for row in csv.reader(csv_file)]
    num_steps = np.round(np.asarray(rows[0], dtype=float))
    times = np.asarray(rows[3], dtype=float) * 1e-9
    num_bisects = np.round(np.asarray(rows[2], dtype=float))

ax_time.plot(num_steps, times)
# ax_time.set_xscale('log')
# ax_time.set_yscale('log')
ax_time.set_title(f'Continuation Times')
ax_time.set_xlabel('Node Buffer Size')
ax_time.set_ylabel('Continuation Time [s]')

ax_num_bisection.plot(num_steps, num_bisects)
# ax_num_bisection.set_xscale('log')
# ax_num_bisection.set_yscale('log')
ax_num_bisection.set_title(f'Number of Bisections')
ax_num_bisection.set_xlabel('Node Buffer Size')
ax_num_bisection.set_ylabel('Number of Bisections')

buffer_sizes, final_nodes = [list(elements) for elements in zip(*sorted(zip(buffer_sizes, final_nodes)))]
ax_final_nodes.plot(buffer_sizes, final_nodes)
# ax_final_nodes.set_xscale('log')
# ax_final_nodes.set_yscale('log')
ax_final_nodes.set_title(f'Final Number of Nodes')
ax_final_nodes.set_xlabel('Node Buffer Size')
ax_final_nodes.set_ylabel('Final Number of Nodes')

# ax_nodes.set_xscale('log')
# ax_nodes.set_yscale('log')
ax_nodes.set_title(f'Num. Nodes over Continuation')
ax_nodes.set_xlabel(r'Terminal Downrange, $\phi_f$ [deg]')
ax_nodes.set_ylabel('Number of Nodes')
fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax_nodes, label='Node Buffer Size')

fig.subplots_adjust(
        top=0.838,
        bottom=0.142,
        left=0.108,
        right=0.936,
        hspace=0.734,
        wspace=0.3
)

plt.show()
