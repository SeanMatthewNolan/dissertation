import numpy as np
import matplotlib.pyplot as plt
import itertools
from matplotlib import rc, rcParams

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

major_num_steps = (4, 4, 4)
minor_num_steps = (3, 3, 3)

lower_bounds = (0, 0, 0)
upper_bounds = (1, 1, 1)

# major_num_steps = (11, 11, 11)
# minor_num_steps = (20, 20, 5)
#
# lower_bounds = (0, 0, 0)
# upper_bounds = (1, 1, 1)

step_list = [lower_bounds]
major_point_cloud = [lower_bounds]
minor_point_cloud = []

for n in range(len(major_num_steps)):
    new_maj_pnts = []
    maj_step_series = np.linspace(lower_bounds[n], upper_bounds[n], major_num_steps[n])
    for maj_point in major_point_cloud:
        maj_pnt_line = []
        for maj_step in maj_step_series:
            new_pnt = list(maj_point)
            new_pnt[n] = maj_step
            maj_pnt_line.append(tuple(new_pnt))

        for maj_pnt_0, maj_pnt_1 in zip(maj_pnt_line[:-1], maj_pnt_line[1:]):
            min_step_series = np.linspace(maj_pnt_0[n], maj_pnt_1[n], minor_num_steps[n]+1, endpoint=False)[1:]
            for min_step in min_step_series:
                new_point = list(maj_pnt_0)
                new_point[n] = min_step
                step_list.append(tuple(new_point))
                minor_point_cloud.append(tuple(new_point))

            step_list.append(maj_pnt_1)

        new_maj_pnts += maj_pnt_line

    major_point_cloud += new_maj_pnts


def plot_pnts(pnts):
    xs, ys, zs = [], [], []
    for pnt in pnts:
        xs.append(pnt[0])
        ys.append(pnt[1])
        zs.append(pnt[2])

    return xs, ys, zs


xx, yy, zz = plot_pnts(step_list)
maj_xx, maj_yy, maj_zz = plot_pnts(major_point_cloud)
min_xx, min_yy, min_zz = plot_pnts(minor_point_cloud)

fig = plt.figure(figsize=(6.5, 3))

ax0 = fig.add_subplot(131, projection='3d')
int_set = [np.linspace(low, upp, (k_maj - 1) + 1)
          for k_maj, k_min, low, upp in zip(major_num_steps, minor_num_steps, lower_bounds, upper_bounds)]
x_full, y_full, z_full = [], [], []
for x, y, z in itertools.product(int_set[0], int_set[1], int_set[2]):
    x_full.append(x)
    y_full.append(y)
    z_full.append(z)
ax0.scatter(x_full, y_full, z_full)
ax0.set_title('Points of Interest')

ax1 = fig.add_subplot(132, projection='3d')
ls_set = [np.linspace(low, upp, (k_maj - 1)*(k_min + 1) + 1)
          for k_maj, k_min, low, upp in zip(major_num_steps, minor_num_steps, lower_bounds, upper_bounds)]
x_sub, y_sub, z_sub = [], [], []
for x, y, z in itertools.product(ls_set[0], ls_set[1], ls_set[2]):
    if not(x in int_set[0] and y in int_set[1] and z in int_set[2]):
        x_sub.append(x)
        y_sub.append(y)
        z_sub.append(z)
ax1.scatter(x_full, y_full, z_full)
ax1.scatter(x_sub, y_sub, z_sub, color='C1')
ax1.set_title('Full Product Space')

ax2 = fig.add_subplot(133, projection='3d')
ax2.scatter(maj_xx, maj_yy, maj_zz, label='Points of Interest')
ax2.scatter(min_xx, min_yy, min_zz, label='Steps for Convergence')
ax2_pos = ax2.get_position().bounds
# ax2.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center')
ax2.set_title('Sparse Product Space')

fig.legend(ncols=2, loc='lower center')

color_var_series = [k for k in range(len(xx))]
# normer = matplotlib.norms.Normalize(vmin=min(color_var_series), vmax=max(color_var_series))
# cm = matplotlib.norms.LinearSegmentedColormap.from_list('mycolors', ['C0', 'C2', 'C1'])
# cm = matplotlib.cm.jet
# sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=normer)
# sm.set_array(color_var_series)
#
# for pnt, color_var in zip(step_list, color_var_series):
#     ax.scatter(pnt[0], pnt[1], pnt[2], color=sm.to_rgba(color_var))

# ax3 = fig.add_subplot(133, projection='3d')
# sc = ax3.scatter(xx, yy, zz, cmap='viridis', c=color_var_series)
# plt.subplots_adjust(bottom=0.2, left=0.1, right=0.9, top=0.9)
# ax3.set_title('Order of Evaulation')
#
# ax3_pos = ax3.get_position().bounds
# cax = plt.axes([ax3_pos[0], ax3_pos[1] - 0.1, ax3_pos[2], 0.025])
# fig.colorbar(sc, cax=cax, format='%g', label='Evaluation Order', orientation='horizontal')

fig.subplots_adjust(
        top=0.95,
        bottom=0.05,
        left=0.031,
        right=0.93,
        hspace=0.2,
        wspace=0.313
)

plt.show()
