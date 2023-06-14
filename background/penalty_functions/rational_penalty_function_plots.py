import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams, rc, gridspec
from matplotlib.colors import LogNorm

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

x_min, x_max = 0, 6


def rational_min(_x):
    return 1 / (_x - x_min)


def rational_max(_x):
    return 1 / (x_max - x)


# size = (3.5, 3.5)
size = (6.5, 4)

x = np.linspace(x_min + 1e-4, x_max - 1e-4, 100)

xticks_val = np.concatenate((np.arange(x_min, x_max, 2), [x_max]))
xticks_lab_min = [r'$C_{i, \text{min}}$'] \
                 + [r'$C_{i, \text{min}} + ' + r'{0:d}$'.format(_val - x_min) for _val in xticks_val[1:]]
xticks_lab_max = [r'$C_{i, \text{max}} - ' + r'{0:d}$'.format(x_max - _val) for _val in xticks_val[:-1]] \
                 + [r'$C_{i, \text{max}}$']

fig = plt.figure(0, figsize=size)
gs = gridspec.GridSpec(2, 2, height_ratios=(1, 0.05))

fig.suptitle('Rational Penalty Functions')

ax_min = fig.add_subplot(gs[0, 0])
ax_max = fig.add_subplot(gs[0, 1])
ax_cb = fig.add_subplot(gs[1, :])

ax_min.set_title('Minimum Value Constraint')
ax_min.set_xticks(xticks_val, xticks_lab_min)
ax_min.set_xlabel(r'$C_i[t, \bm{x}{(t)}, \bm{u}{(t)}]$')
ax_min.set_ylabel(r'Penalty Cost')

ax_max.set_title('Maximum Value Constraint')
ax_max.set_xticks(xticks_val, xticks_lab_max)
ax_max.set_xlabel(r'$C_i[t, \bm{x}{(t)}, \bm{u}{(t)}]$')
ax_max.set_ylabel(r'Penalty Cost')

cmap = plt.colormaps['viridis_r']

epsilons = np.geomspace(1e-2, 1, 10)
norm = LogNorm(epsilons[0], epsilons[-1])
for eps in epsilons:
    color = cmap(norm(eps))
    ax_min.plot(x, eps * rational_min(x), color=color)
    ax_max.plot(x, eps * rational_max(x), color=color)

ax_min.set_xlim(x_min, x_max)
ax_min.set_ylim(0, 10)

ax_max.set_xlim(x_min, x_max)
ax_max.set_ylim(0, 10)

fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax_cb, orientation='horizontal',
             label=r'Smoothing Constant, $\epsilon$')

fig.subplots_adjust(
        top=0.838,
        bottom=0.137,
        left=0.087,
        right=0.948,
        hspace=0.516,
        wspace=0.365
)

plt.show()
