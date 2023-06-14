import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

cmap = plt.colormaps['viridis']

sols = gp.load_sol_set('sol_set_alg.bin')

fig = plt.figure(figsize=(6.5, 3.5))
title = fig.suptitle('Adjoint Variables and Feasibility Boundary: Brachistochrone')

gs = gridspec.GridSpec(2, 2, height_ratios=(1, 0.1))

ax_path = fig.add_subplot(gs[0, 0])
ax_norm = fig.add_subplot(gs[0, 1])
# ax_adjoint = fig.add_subplot(gs[0, 2])
ax_cb_v = fig.add_subplot(gs[1, :])

yfs = []
nus = []
adj_norm = []
for sol in sols:
    yfs.append(sol.x[1, -1])
    nus.append(np.concatenate((sol.nu0, sol.nuf)))
    adj_norm.append(np.linalg.norm(np.concatenate((sol.nu0, sol.nuf))))

yfs = np.array(yfs)
nus = np.array(nus)

# ax_adjoint.plot(yfs, nus)
# ax_adjoint.set_title('Static Adjoint Variables')
# ax_adjoint.set_xlabel('Terminal Vertical Position, $y_f$ [ft]')
# ax_adjoint.set_ylabel(r'Static Adjoints Norm, $|| \boldsymbol{\nu} ||$')

ax_norm.plot(yfs, np.array(adj_norm))
ax_norm.set_title('Norm of Static Adjoint Variables')
ax_norm.set_xlabel('Terminal Vertical Position, $y_f$ [ft]')
ax_norm.set_ylabel(r'Static Adjoints Norm, $|| \boldsymbol{\nu} ||$')

norm = Normalize(min(yfs), max(yfs))

for sol in sols[::]:
    ax_path.plot(sol.x[0, :], sol.x[1, :], color=cmap(norm(sol.x[1, -1])), linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax_cb_v, label='Terminal Vertical Position [ft]', orientation='horizontal', aspect=50)

ax_path.set_title('Paths')
ax_path.set_xlabel('$x$ [ft]')
ax_path.set_ylabel('$y$ [ft]')

fig.subplots_adjust(
        top=0.815,
        bottom=0.161,
        left=0.122,
        right=0.977,
        hspace=0.662,
        wspace=0.24
)

plt.show()
