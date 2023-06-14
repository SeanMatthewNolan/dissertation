import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc
from matplotlib.colors import Normalize

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

# max_phi = gp.load_sol_set('max_phi.bin')
max_v = gp.load_sol_set('max_v.bin')
max_v_r = gp.load_sol_set('max_v_reverse.bin')

max_v.solutions = max_v_r.solutions[::-1] + max_v.solutions

fig = plt.figure(figsize=(3.5, 2.9))
title = fig.suptitle('Feasibility Boundary:\nHGV Maximum Terminal Velocity')

gs = gridspec.GridSpec(2, 1)

ax_adjoint = fig.add_subplot(gs[1, 0])
ax_max_v = fig.add_subplot(gs[0, 0])

max_v_phi = []
max_v_v = []
nus = []
adj_norm = []
for sol in max_v:
    max_v_phi.append(sol.x[1, -1])
    max_v_v.append(sol.x[2, -1])
    nus.append(np.concatenate((sol.nu0, sol.nuf)))
    adj_norm.append(np.linalg.norm(np.concatenate((sol.nu0, sol.nuf))))

nus = np.array(nus)

cmap = plt.colormaps['viridis']
cmap_v = plt.colormaps['viridis_r']
max_phif = np.rad2deg(max_v[-1].x[1, -1])
min_phif = np.rad2deg(max_v[0].x[1, -1])

ax_adjoint.plot(np.rad2deg(max_v_phi), np.array(adj_norm))
ax_adjoint.ticklabel_format(scilimits=(-5, 3))

# ax_adjoint.set_title('Norm of Static Adjoint Variables')
# ax_adjoint.set_yscale('log')
ax_adjoint.set_xlabel('Terminal Downrange Angle, $\phi_f$ [deg]')
ax_adjoint.set_ylabel(r'$\left\lVert\bm{\nu}\right\lVert$')

for sol in max_v[::]:
    ax_max_v.plot(np.rad2deg(sol.x[1, :]), sol.x[0, :] / 1e5,
                  color=cmap_v((np.rad2deg(sol.x[1, -1]) - min_phif) / (max_phif - min_phif)), linewidth=1)

fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(min_phif, max_phif), cmap=cmap_v), ax=ax_max_v,
                  label=r'$\phi_f$ [deg]', orientation='vertical')

# ax_max_v.set_title('Max Terminal Velocity')
ax_max_v.set_xlabel('Downrange Angle, $\phi$ [deg]')
ax_max_v.set_ylabel('Altitude,\n' + r'$h$ [$\times 10^5$ ft]')

fig.subplots_adjust(
        top=0.779,
        bottom=0.196,
        left=0.19,
        right=0.93,
        hspace=1.0,
        wspace=0.2
)

plt.show()
