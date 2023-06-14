import pickle
from math import log

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

cmap = plt.colormaps['viridis']

fig = plt.figure(figsize=(6.5, 7.5))
fig.suptitle('Performance of New Guess Generation Method:\nCart Navigation Problem Beacon Position')

gs = gridspec.GridSpec(4, 2)

ax_path_05 = fig.add_subplot(gs[0, 0])
ax_path_10 = fig.add_subplot(gs[1, 0])
ax_path_25 = fig.add_subplot(gs[2, 0])
ax_path_50 = fig.add_subplot(gs[3, 0])

ax_term_05 = fig.add_subplot(gs[0, 1])
ax_term_10 = fig.add_subplot(gs[1, 1])
ax_term_25 = fig.add_subplot(gs[2, 1])
ax_term_50 = fig.add_subplot(gs[3, 1])

axs = [ax_path_05, ax_path_10, ax_path_25, ax_path_50, ax_term_05, ax_term_10, ax_term_25, ax_term_50]
files = ['nav_path_t5.data', 'nav_path_t10.data', 'nav_path_t25.data', 'nav_path_t50.data',
         'nav_final_t5.data', 'nav_final_t10.data', 'nav_final_t25.data', 'nav_final_t50.data']
titles = [r'Path Cost: $t_{\text{span}}$ = ' + f'{t_span:d} [sec]' for t_span in [5, 10, 25, 50]] \
         + [r'Terminal Cost: $t_{\text{span}}$ = ' + f'{t_span:d} [sec]' for t_span in [5, 10, 25, 50]]

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


for filename, ax, title in zip(files, axs, titles):

    with open(filename, 'rb') as file:
        data = pickle.load(file)

    n_conv = 0
    for data_dict in data:
        guess = data_dict['guess']
        converged = data_dict['converged']
        xb, yb = data_dict['xb'], data_dict['yb']

        marker, marker_color, label = _gen_markers(converged)
        if converged:
            n_conv += 1

        ax.plot(xb, yb, marker=marker, markersize=markersize, markeredgewidth=markeredgewidth,
                color=marker_color, label=label, linestyle='')

    ax.set_title(title)
    ax.set_xlabel('$x_b$ [ft]')
    ax.set_ylabel('$y_b$ [ft]')

    print(f'{filename}: {n_conv / len(data) * 100:4.3f} %')

handles, labels = axs[2].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(), markerscale=4 / markersize, loc='lower center', ncols=2)

fig.subplots_adjust(
        top=0.89,
        bottom=0.105,
        left=0.1,
        right=0.977,
        hspace=0.814,
        wspace=0.257
)

plt.show()
