import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec, rc

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

alg_set = gp.load_sol_set('downrange_alg.data')
dif_set = gp.load_sol_set('downrange_dif.data')

alg_sol = alg_set[-1]
dif_sol = dif_set[-1]

fig = plt.figure(figsize=(6.5, 3))
title = fig.suptitle('Space Shuttle Maximum Downrange')

gs = gridspec.GridSpec(1, 3)

ax_path = fig.add_subplot(gs[0, 0])
ax_hv = fig.add_subplot(gs[0, 1])
ax_alpha = fig.add_subplot(gs[0, 2])

ax_path.plot(alg_sol.x[1, :] * 180 / np.pi, alg_sol.x[0, :] / 1e5, '-',
             label='Algebraic Control Law')
ax_path.plot(dif_sol.x[1, :] * 180 / np.pi, dif_sol.x[0, :] / 1e5, '--',
             label=r'Differential Control Law')
ax_path.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_path.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

fig.legend(loc='lower center', ncols=2)

ax_hv.plot(alg_sol.x[2, :] / 1e3, alg_sol.x[0, :] / 1e5, '-',  label='Algebraic Control Law')
ax_hv.plot(dif_sol.x[2, :] / 1e3, dif_sol.x[0, :] / 1e5, '--', label=r'Differential Control Law')
ax_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
ax_hv.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')


ax_alpha.plot(alg_sol.t, alg_sol.u[0, :] * 180 / np.pi, '-',  label='Algebraic Control Law')
ax_alpha.plot(dif_sol.t, dif_sol.u[0, :] * 180 / np.pi, '--', label=r'Differential Control Law')
ax_alpha.set_xlabel(r'Time $t$ [s]')
ax_alpha.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

fig.subplots_adjust(
        top=0.85,
        bottom=0.3,
        left=0.1,
        right=0.975,
        hspace=0.2,
        wspace=0.6
)

plt.show()
