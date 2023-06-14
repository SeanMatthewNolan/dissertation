import pickle
from math import log

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

# cmap_alpha = plt.colormaps['viridis']
# cmap_alt = plt.colormaps['cividis']
cmap_down = plt.colormaps['plasma']

cmap_alt = plt.colormaps['viridis']
# cmap_down = plt.colormaps['viridis']
# cmap_vel = plt.colormaps['viridis']
# cmap_fpa = plt.colormaps['viridis']

file_name = 'giuseppe_guess_test.pickle'
# file_name = 'giuseppe_guess_test_simpson.pickle'
# file_name = 'giuseppe_guess_test_midpoint.pickle'

with open(file_name, 'rb') as file:
    data = pickle.load(file)

fig = plt.figure(figsize=(6, 5.5))
fig.suptitle('Performance of New Guess Generation Method:\nSpace Shuttle Maximum Crossrange Problem')

gs0 = gridspec.GridSpec(2, 1, height_ratios=[0.05, 2.05], top=1-0.1, hspace=0.2)
gs_ldg = gs0[0, :].subgridspec(1, 2, wspace=1)
gs_plt = gs0[1:, :].subgridspec(3, 2, height_ratios=[1, 1, 0.1], hspace=0.75, wspace=0.5)

ax_conv = fig.add_subplot(gs_plt[0, 0])
ax_path = fig.add_subplot(gs_plt[0, 1])

ax_alt = fig.add_subplot(gs_plt[1, 0])
cs_ax_alt = fig.add_subplot(gs_plt[2, 0])

ax_down = fig.add_subplot(gs_plt[1, 1])
cs_ax_down = fig.add_subplot(gs_plt[2, 1])

alphas = data[2]
max_alpha = max(alphas)
min_alpha = min(alphas)

alts = [guess.x[0, -1] for guess in data[0]]
downs = [guess.x[1, -1] for guess in data[0]]
vels = [guess.x[3, -1] for guess in data[0]]
fpas = [guess.x[4, -1] for guess in data[0]]

# min_alt = min(alts)
max_alt = max(alts)

min_alt = 100_000
# max_alt = 300_000

min_down = min(downs)
max_down = max(downs)

min_vel = min(vels)
max_vel = max(vels)

min_fpa = -15 / 180 * 3.14159
max_fpa = max(fpas)


def scale_alpha(_var):
    return (_var - min_alpha) / (max_alpha - min_alpha)


def scale_alt(_var):
    return (_var - min_alt) / (max_alt - min_alt)


def scale_down(_var):
    return (log(_var) - log(min_down)) / (log(max_down) - log(min_down))


def scale_vel(_var):
    return (_var - min_vel) / (max_vel - min_vel)


def scale_fpa(_var):
    return (_var - min_fpa) / (max_fpa - min_fpa)


linewidth = 0.5
markersize = 2
markeredgewidth = 0.25
for guess, t_span, alpha, converged in zip(*data):
    alt = guess.x[0, -1]
    down = guess.x[1, -1]
    vel = guess.x[3, -1]
    fpa = guess.x[4, -1]

    if converged:
        marker = 'o'
        marker_color = 'C0'
        label = 'Converged'
    else:
        marker = 'x'
        marker_color = 'C1'
        label = 'Not Converged'

    # color = cmap_alpha(scale_alpha(alpha))
    # ax_path.plot(guess.x[1, :] * 180 / 3.14159, guess.x[0, :] / 1000,
    #              color=cmap_alpha(scale_alpha(alpha)), linewidth=linewidth)
    ax_path.plot(guess.x[1, :] * 180 / 3.14159, guess.x[0, :] / 1000,
                 color=marker_color, label=label, linewidth=linewidth)

    ax_alt.plot(alpha, t_span, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
                color=cmap_alt(scale_alt(alt)), label=label, linestyle='')
    ax_down.plot(alpha, t_span, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
                 color=cmap_down(scale_down(down)), label=label, linestyle='')
    # ax_vel.plot(alpha, t_span, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
    #             color=cmap_vel(scale_vel(vel)), label=label, linestyle='')
    # ax_fpa.plot(alpha, t_span, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
    #             color=cmap_fpa(scale_fpa(fpa)), label=label, linestyle='')
    ax_conv.plot(alpha, t_span, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
                 color=marker_color, label=label, linestyle='')


print(f'Converged {sum(data[3])/len(data[3])*100:3.1f}%')


ax_path.set_xlabel(r'Downrange Angle, $\phi$ [deg]')
# ax_path.set_ylabel(r'Crossrange Angle, $\theta$ [deg]')
# ax_path.set_zlabel(r'Altitude, $h$ [ft]')
ax_path.set_ylabel(r'Altitude, $h$ [$\times10^3$ ft]')
# x_lim = max(ax_path.get_xlim())
# ax_path.set_ylim(-x_lim/2, x_lim/2)

ax_alt.set_yscale('log')
ax_alt.set_xlabel(r'AoA, $\alpha$ [deg]')
ax_alt.set_ylabel(r'Final Time, $t_f$ [s]')

ax_down.set_yscale('log')
ax_down.set_xlabel(r'AoA, $\alpha$ [deg]')
ax_down.set_ylabel(r'Final Time, $t_f$ [s]')

# ax_vel.set_yscale('log')
# ax_vel.set_xlabel(r'AoA, $\alpha$ [deg]')
# ax_vel.set_ylabel(r'Final Time, $t_f$ [s]')

# ax_fpa.set_yscale('log')
# ax_fpa.set_xlabel(r'AoA, $\alpha$ [deg]')
# ax_fpa.set_ylabel(r'Final Time, $t_f$ [s]')

ax_conv.set_yscale('log')
# ax_conv.set_ylim((None, 1000))
ax_conv.set_xlabel(r'AoA, $\alpha$ [deg]')
ax_conv.set_ylabel(r'Final Time, $t_f$ [s]')

handles, labels = ax_conv.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(), markerscale=4 / markersize, loc='center',
           bbox_to_anchor=gs_ldg[0, 0].get_position(fig), ncols=2, fontsize=8, frameon=False)

handles, labels = ax_path.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(),
           bbox_to_anchor=gs_ldg[0, 1].get_position(fig), loc='center', fontsize=8, ncols=2, frameon=False, )

orientation = 'horizontal'
# fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_alpha, max_alpha), cmap=cmap_alpha),
#              orientation=orientation, aspect=50,
#              cax=cs_ax_path, label=r'AoA in Guess, $\alpha$ [deg]')

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_alt / 1000, max_alt / 1000), cmap=cmap_alt),
             orientation=orientation, aspect=100, extend='min',
             cax=cs_ax_alt, label=r'Term. Altitude, $h_f$ [$\times10^3$ ft]')

fig.colorbar(plt.cm.ScalarMappable(norm=LogNorm(min_down * 180 / 3.14159, max_down * 180 / 3.14159), cmap=cmap_down),
             orientation=orientation, aspect=100,
             cax=cs_ax_down, label=r'Term. Downrange, $\phi_f$ [deg]')

# fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_vel / 1000, max_vel / 1000), cmap=cmap_vel),
#              orientation=orientation, aspect=50,
#              cax=cs_ax_vel, label=r'Term. Velocity, $v_f$ [$\times10^3$ ft/s]')

# fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_fpa * 180 / 3.14159, max_fpa * 180 / 3.14159), cmap=cmap_fpa),
#              orientation=orientation, aspect=50, extend='min',
#              cax=cs_ax_fpa, label=r'Term. FPA, $\gamma_f$ [deg]')

gs0.tight_layout(fig)

plt.show()
