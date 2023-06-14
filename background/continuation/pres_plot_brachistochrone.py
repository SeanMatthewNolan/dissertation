import matplotlib.pyplot as plt
import giuseppe as gp
from matplotlib import rcParams, gridspec, rc

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

cmap = plt.colormaps['viridis']

sol_set_alg = gp.data_classes.load_sol_set('sol_set_alg.bin')

fig_cont = plt.figure(figsize=(2.75, 2.75))
ax_cont = fig_cont.add_subplot()

xf, yf = [], []
for sol in sol_set_alg:
    ax_cont.plot(sol.x[0, :], sol.x[1, :], color='C0')
    xf.append(sol.x[0, -1])
    yf.append(sol.x[1, -1])
ax_cont.plot(xf, yf, '--', marker='*', label='', color='C1')

offset, mult = 0.25, 0.95
for xfl, xfr, yfl, yfr in zip(xf[:-1], xf[1:], yf[:-1], yf[1:]):
    ax_cont.arrow(xfl + offset, yfl + offset, mult * (xfr - xfl), mult * (yfr - yfl),
                  width=0.1, length_includes_head=True, linewidth=0, color='C1')

ax_cont.set_title('Brachistochrone Problem')
ax_cont.set_aspect('equal', 'box')
ax_cont.set_xlabel(r'$x$ [ft]')
ax_cont.set_ylabel(r'$y$ [ft]')

fig_cont.tight_layout()

plt.show()
