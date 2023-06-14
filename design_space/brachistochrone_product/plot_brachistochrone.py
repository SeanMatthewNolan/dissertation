import matplotlib.pyplot as plt
import numpy as np

import giuseppe as gp
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

cmap = plt.colormaps['viridis']

# sol_set_alg = gp.data_classes.load_sol_set('sol_set_alg.bin')
# sol_set_diff = gp.data_classes.load_sol_set('sol_set_diff.bin')
sol_set = gp.data_classes.load_sol_set('sol_set_product.bin')

fig_cont = plt.figure(figsize=(6.5, 3))
ax_cont = fig_cont.add_subplot()

xf, yf = [], []
for sol in sol_set:
    ax_cont.plot(sol.x[0, :], sol.x[1, :], color='C0')
    xf.append(sol.x[0, -1])
    yf.append(sol.x[1, -1])
ax_cont.plot(xf, yf, '*', label='', color='C1')

offset, mult = 0.25, 0.95

x_span = np.linspace(10, 100, 10)
for xfl, xfr in zip(x_span[:-1], x_span[1:]):
    ax_cont.arrow(xfl, 0, mult * (xfr - xfl), 0,
                  width=1, length_includes_head=True, linewidth=0, color='C1')

y_span = -np.linspace(10, 100, 10)
for yfl, yfr in zip(y_span[:-1], y_span[1:]):
    ax_cont.arrow(0, yfl, 0, mult * (yfr - yfl),
                  width=1, length_includes_head=True, linewidth=0, color='C1', zorder=100)

ax_cont.set_title('Product Space Example:\nBrachistochrone Problem')
ax_cont.set_aspect('equal', 'box')
ax_cont.set_xlabel(r'$x$ [ft]')
ax_cont.set_ylabel(r'$y$ [ft]')


fig_cont.subplots_adjust(
        top=0.852,
        bottom=0.161,
        left=0.023,
        right=0.977,
        hspace=0.2,
        wspace=0.2
)

plt.show()
