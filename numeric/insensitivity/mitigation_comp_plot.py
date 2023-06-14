import giuseppe as gp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, gridspec

rcParams['font.family'] = 'serif'
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath,bm}'

alg_set = gp.load_sol_set('downrange_alg_10.data')
dif_set = gp.load_sol_set('downrange_dif_mit.data')
# dif_set = gp.load_sol_set('downrange_dif_mit.data')

alg_sol = alg_set[-1]
dif_sol = dif_set[-1]

fig = plt.figure(figsize=(6.5, 4.25))
title = fig.suptitle('Mitigation of Control Insensitivity')

gs = gridspec.GridSpec(2, 2)

ax_path = fig.add_subplot(gs[0, 0])
ax_hv = fig.add_subplot(gs[0, 1])
ax_alpha = fig.add_subplot(gs[1, 0])
ax_huu = fig.add_subplot(gs[1, 1])

ax_path.plot(dif_sol.x[1, :] * 180 / np.pi, dif_sol.x[0, :] / 1e5, '-',
             label=r'Differential Control with Mitigation')
ax_path.plot(alg_sol.x[1, :] * 180 / np.pi, alg_sol.x[0, :] / 1e5, '--',
             label='Algebraic Control')
ax_path.set_xlabel(r'Downrange, $\phi$ [deg]')
ax_path.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')

fig.legend(loc='lower center', ncols=2)

ax_hv.plot(dif_sol.x[2, :] / 1e3, dif_sol.x[0, :] / 1e5, '-', label=r'Differential Control with Mitigation')
ax_hv.plot(alg_sol.x[2, :] / 1e3, alg_sol.x[0, :] / 1e5, '--',  label='Algebraic Control')
# ax_hv.set_yscale('log')
ax_hv.set_ylim(0, 5)
ax_hv.set_xlabel(r'Velocity, $v$ [$\times10^3$ ft/s]')
ax_hv.set_ylabel(r'Altitude, $h$ [$\times10^5$ ft]')


ax_alpha.plot(dif_sol.t, dif_sol.u[0, :] * 180 / np.pi, '-', label=r'Differential Control with Mitigation')
ax_alpha.plot(alg_sol.t, alg_sol.u[0, :] * 180 / np.pi, '--',  label='Algebraic Control')
ax_alpha.set_xlabel(r'Time $t$ [s]')
ax_alpha.set_ylabel(r'Angle-of-Attack, $\alpha$ [deg]')

ax_huu.plot(dif_sol.t, dif_sol.aux['H_uu'][:, 0, 0], '-', label=r'Differential Control with Mitigation')
ax_huu.plot(alg_sol.t, alg_sol.aux['H_uu'][:, 0, 0], '--',  label='Algebraic Control')
ax_huu.set_xlabel(r'Time $t$ [s]')
ax_huu.set_ylabel(r'$H_{\bm{uu}}$ [1/rad/s]')
ax_huu.set_yscale('log')

fig.subplots_adjust(
        top=0.895,
        bottom=0.188,
        left=0.091,
        right=0.977,
        hspace=0.43,
        wspace=0.311
)

print(f'{np.rad2deg(alg_sol.x[1, -1] - dif_sol.x[1, -1])}')
print(f'{(alg_sol.x[1, -1] - dif_sol.x[1, -1]) / alg_sol.x[1, -1] * 100}')

plt.show()
