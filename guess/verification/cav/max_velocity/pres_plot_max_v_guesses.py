import pickle
from math import log

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

cmap = plt.colormaps['viridis']

cav_l_file_name = 'cav_l_data.data'

with open(cav_l_file_name, 'rb') as file:
    cav_l_data = pickle.load(file)

cav_h_file_name = 'cav_h_data.data'

with open(cav_h_file_name, 'rb') as file:
    cav_h_data = pickle.load(file)

fig = plt.figure(figsize=(3.5, 2.9))

gs = gridspec.GridSpec(2, 1)

ax_cav_l_conv = fig.add_subplot(gs[0, 0])
ax_cav_l_traj = fig.add_subplot(gs[1, 0])
# ax_cav_l_hv = fig.add_subplot(gs[2, 0])

# ax_cav_h_conv = fig.add_subplot(gs[0, 1])
# ax_cav_h_traj = fig.add_subplot(gs[1, 1])
# ax_cav_h_hv = fig.add_subplot(gs[2, 1])

linewidth = 0.5
markersize = 2
markeredgewidth = 0.25


def _gen_markers(_converged):
    if converged:
        _marker = 'o'
        _marker_color = 'C0'
        _label = 'Converged'
    else:
        _marker = 'x'
        _marker_color = 'C1'
        _label = 'Not Converged'

    return _marker, _marker_color, _label


n_conv = 0
for data_dict in cav_l_data:
    guess = data_dict['guess']
    converged = data_dict['converged']
    h0, v0 = data_dict['h0'], data_dict['v0']

    marker, marker_color, label = _gen_markers(converged)
    if converged:
        n_conv += 1

    ax_cav_l_conv.plot(v0 / 1000, h0 / 1000, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
                       color=marker_color, label=label, linestyle='')
    ax_cav_l_traj.plot(guess.x[1, :] * 180 / 3.14159, guess.x[0, :] / 1000,
                       color=marker_color, label=label, linewidth=linewidth)
    # ax_cav_l_hv.plot(guess.x[3, :] / 1000, guess.x[0, :] / 1000,
    #                  color=marker_color, label=label, linewidth=linewidth)


print(f'CAV-L Converge: {n_conv / len(cav_l_data) * 100:4.3f} %')


# ax_cav_l_conv.set_title('CAV-L: Guess Constants')
ax_cav_l_conv.set_xlabel(r'Init.\ Velocity, $v_0$ [$\times10^3$ ft/s]')
ax_cav_l_conv.set_ylabel(r'Init.\ Altitude, ' '\n' r'$h_0$ [$\times10^3$ ft]')

# ax_cav_l_traj.set_title('CAV-L: Guess Paths')
ax_cav_l_traj.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_cav_l_traj.set_ylabel(r'Altitude, $h$ ' '\n' r'[$\times10^3$ ft]')

# ax_cav_l_hv.set_title('CAV-L: Guess h-v')
# ax_cav_l_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
# ax_cav_l_hv.set_ylabel(r'Altitude, $h$ [$\times10^3$ ft]')

n_conv = 0
for data_dict in cav_h_data:
    guess = data_dict['guess']
    converged = data_dict['converged']
    h0, v0 = data_dict['h0'], data_dict['v0']

    marker, marker_color, label = _gen_markers(converged)
    if converged:
        n_conv += 1

    # ax_cav_h_conv.plot(v0 / 1000, h0 / 1000, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
    #                    color=marker_color, label=label, linestyle='')
    # ax_cav_h_traj.plot(guess.x[1, :] * 180 / 3.14159, guess.x[0, :] / 1000,
    #                    color=marker_color, label=label, linewidth=linewidth)
    # ax_cav_h_hv.plot(guess.x[3, :] / 1000, guess.x[0, :] / 1000,
    #                  color=marker_color, label=label, linewidth=linewidth)

print(f'CAV-H Converge: {n_conv / len(cav_h_data) * 100:4.3f} %')

# ax_cav_h_conv.set_title('CAV-H: Guess Constants')
# ax_cav_h_conv.set_xlabel(r'Init.\ Velocity, $v_0$ [$\times10^3$ ft/s]')
# ax_cav_h_conv.set_ylabel(r'Init.\ Altitude, $h_0$ [$\times10^3$ ft]')
#
# ax_cav_h_traj.set_title('CAV-H: Guess Paths')
# ax_cav_h_traj.set_xlabel(r'Downrange, $\phi$ [deg]')
# ax_cav_h_traj.set_ylabel(r'Altitude, $h$ [$\times10^3$ ft]')
#
# ax_cav_h_hv.set_title('CAV-H: Guess h-v')
# ax_cav_h_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
# ax_cav_h_hv.set_ylabel(r'Altitude, $h$ [$\times10^3$ ft]')

handles, labels = ax_cav_l_conv.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(), markerscale=4 / markersize, loc='lower center', ncols=2)

fig.subplots_adjust(
        top=0.948,
        bottom=0.28,
        left=0.229,
        right=0.957,
        hspace=0.711,
        wspace=0.261
)

plt.show()
