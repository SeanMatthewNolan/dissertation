import pickle
from math import log

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams, gridspec, rc
import numpy as np

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

cmap = plt.colormaps['viridis']
# cmap = plt.colormaps['magma']

# file_name = 'beluga_guess_test_10s.pickle'
file_name = 'beluga_guess_test_100s.pickle'
# file_name = 'beluga_guess_test_200s.pickle'

with open(file_name, 'rb') as file:
    data = pickle.load(file)

fig = plt.figure(figsize=(3.5, 2.9))
fig.suptitle('Performance of Former\nGuess Generation Method')


# gs = gridspec.GridSpec(6, 3, width_ratios=[1, 1, 0.025], height_ratios=[2/3, 2/3, 2/3, 2/3, 2/3, 2/3])
# ax_error = fig.add_subplot(gs[0:3, 0])
# ax_ratio = fig.add_subplot(gs[3:6, 0])
# ax_path = fig.add_subplot(gs[0:2, 1])
# ax_alpha = fig.add_subplot(gs[2:4, 1])
# ax_beta = fig.add_subplot(gs[4:6, 1])
# ax_color = fig.add_subplot(gs[:, 2])

gs = gridspec.GridSpec(2, 2, height_ratios=[1, 0.05])
# ax_error = fig.add_subplot(gs[0, 0])
ax_ratio = fig.add_subplot(gs[0, 0])
ax_path = fig.add_subplot(gs[0, 1])
# ax_alpha = fig.add_subplot(gs[1, 2:4])
# ax_beta = fig.add_subplot(gs[1, 4:6])
ax_color = fig.add_subplot(gs[1, :])

variations = data[1]
max_var = max(variations)
min_var = min(variations)

converged_dict = {}


def scale_var(_var):
    return (log(_var) - log(min_var)) / (log(max_var) - log(min_var))


linewidth = 0.5
markersize = 2
markeredgewidth = 0.25
for guess, variation, bc_res, converged in zip(*data):
    color = cmap(scale_var(variation))
    # ax_path.plot(guess.x[1, :] * 180 / 3.14159,  guess.x[2, :] * 180 / 3.14159, guess.x[0, :])
    ax_path.plot(guess.x[1, :] * 180 / 3.14159, guess.x[0, :] / 1000, color=color, linewidth=linewidth)
    # ax_alpha.plot(guess.t, guess.u[0, :] * 180 / 3.14159, color=color, linewidth=linewidth)
    # ax_beta.plot(guess.t, guess.u[1, :] * 180 / 3.14159, color=color, linewidth=linewidth)

    if converged:
        marker = 'o'
        marker_color = 'C0'
        label = 'Converged'
    else:
        marker = 'x'
        marker_color = 'C1'
        label = 'Not Converged'

    # ax_error.plot(variation, bc_res, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
    #               color=marker_color, label=label, linestyle='')

    try:
        converged_dict[variation].append(converged)
    except KeyError:
        converged_dict[variation] = [converged]

conv_variations = []
conv_ratios = []
for key, conv_list in converged_dict.items():
    conv_variations.append(key)
    conv_ratios.append(sum(conv_list) / len(conv_list))


conv_variations = np.array(conv_variations)
conv_ratios = np.array(conv_ratios)

ax_ratio.set_xscale('log')
# ax_ratio.plot(conv_variations, conv_ratios * 100)
ax_ratio.bar(conv_variations, conv_ratios * 100, conv_variations * 0.75 * (log(max_var) - log(min_var)) / len(conv_ratios),
             color=[cmap(scale_var(variation)) for variation in conv_variations])
ax_ratio.set_xlabel(r'Perturbation Magnitude, $m$')
ax_ratio.set_ylabel(r'\% Converged')

ax_path.set_xlabel(r'Downrange, $\phi$ [deg]')
# ax_path.set_ylabel(r'Crossrange Angle, $\theta$ [deg]')
# ax_path.set_zlabel(r'Altitude, $h$ [ft]')
ax_path.set_ylabel(r'Altitude, $h$ [$\times10^3$ ft]')
# x_lim = max(ax_path.get_xlim())
# ax_path.set_ylim(-x_lim/2, x_lim/2)

# ax_alpha.set_xlabel(r'Time, $t$ [s]')
# ax_alpha.set_ylabel(r'AoA, $\alpha$ [deg]')
# ax_alpha.set_ylim(-180, 180)
#
# ax_beta.set_xlabel(r'Time, $t$ [s]')
# ax_beta.set_ylabel(r'Bank, $\beta$ [deg]')
# ax_beta.set_ylim(-90, 90)

# ax_error.scatter(data[1], data[2])
# ax_error.set_xscale('log')
# ax_error.set_yscale('log')
# ax_error.set_xlabel(r'Perturbation Magnitude, $m$')
# ax_error.set_ylabel(r'BC Residual')
# handles, labels = ax_error.get_legend_handles_labels()
# by_label = dict(zip(labels, handles))
# fig.legend(by_label.values(), by_label.keys(), markerscale=2, fontsize=8, loc='upper left')

fig.colorbar(plt.cm.ScalarMappable(norm=LogNorm(min_var, max_var), cmap=cmap),
             cax=ax_color, orientation='horizontal',
             label=r'Initial Costates Perturbation Magnitude')

# fig.tight_layout()

plt.show()
